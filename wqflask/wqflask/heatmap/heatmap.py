from __future__ import absolute_import, print_function, division

import sys
sys.path.append(".")

import gc
import string
import cPickle
import os
import datetime
import time
import pp
import math
import collections
import resource

import scipy
import numpy as np
from scipy import linalg

from pprint import pformat as pf

from htmlgen import HTMLgen2 as HT
import reaper

from base.trait import GeneralTrait
from base import data_set
from base import species
from base import webqtlConfig
from utility import webqtlUtil
from wqflask.my_pylmm.data import prep_data
from wqflask.my_pylmm.pyLMM import lmm
from wqflask.my_pylmm.pyLMM import input
from utility import helper_functions
from utility import Plot, Bunch
from utility import temp_data

from MySQLdb import escape_string as escape

import cPickle as pickle
import simplejson as json

from pprint import pformat as pf

from redis import Redis
Redis = Redis()

from flask import Flask, g

class Heatmap(object):

    def __init__(self, start_vars, temp_uuid):
    
        trait_db_list = [trait.strip() for trait in start_vars['trait_list'].split(',')]
        
        helper_functions.get_trait_db_obs(self, trait_db_list)
        
        self.dataset = self.trait_list[0][1]
        
        self.json_data = {} #The dictionary that will be used to create the json object that contains all the data needed to create the figure
        
        self.all_sample_list = []
        self.traits = []
        for trait_db in self.trait_list:
            this_trait = trait_db[0]
            self.traits.append(this_trait.name)
            this_sample_data = this_trait.data
            
            for sample in this_sample_data:
                if sample not in self.all_sample_list:
                    self.all_sample_list.append(sample)
                    
        self.sample_data = []
        for trait_db in self.trait_list:
            this_trait = trait_db[0]
            this_sample_data = this_trait.data
            
            #self.sample_data[this_trait.name] = []
            this_trait_vals = []
            for sample in self.all_sample_list:
                if sample in this_sample_data:
                    this_trait_vals.append(this_sample_data[sample].value)
                    #self.sample_data[this_trait.name].append(this_sample_data[sample].value)
                else:
                    this_trait_vals.append('')
                    #self.sample_data[this_trait.name].append('')
            self.sample_data.append(this_trait_vals)

        self.trait_results = {}
        for trait_db in self.trait_list:
            this_trait = trait_db[0]
            #this_db = trait_db[1]
            self.dataset.group.get_markers()
            
            this_db_samples = self.dataset.group.samplelist
            this_sample_data = this_trait.data
            #print("this_sample_data", this_sample_data)
            this_trait_vals = []
            for index, sample in enumerate(this_db_samples):
                if sample in this_sample_data:
                    sample_value = this_sample_data[sample].value
                    this_trait_vals.append(sample_value)
                else:
                    this_trait_vals.append("x")
                    
            pheno_vector = np.array([val == "x" and np.nan or float(val) for val in this_trait_vals])
            
            key = "pylmm:input:" + str(temp_uuid)
            #print("key is:", pf(key))
            
            genotype_data = [marker['genotypes'] for marker in self.dataset.group.markers.markers]
            
            no_val_samples = self.identify_empty_samples(this_trait_vals)
            trimmed_genotype_data = self.trim_genotypes(genotype_data, no_val_samples)
            
            genotype_matrix = np.array(trimmed_genotype_data).T
            
            #print("genotype_matrix:", str(genotype_matrix.tolist()))
            #print("pheno_vector:", str(pheno_vector.tolist()))
            
            params = dict(pheno_vector = pheno_vector.tolist(),
                        genotype_matrix = genotype_matrix.tolist(),
                        restricted_max_likelihood = True,
                        refit = False,
                        temp_uuid = str(temp_uuid),
                        
                        # meta data
                        timestamp = datetime.datetime.now().isoformat(),
                        )
            
            json_params = json.dumps(params)
            #print("json_params:", json_params)
            Redis.set(key, json_params)
            Redis.expire(key, 60*60)
            print("before printing command")
            
            command = 'python /home/zas1024/gene/wqflask/wqflask/my_pylmm/pyLMM/lmm.py --key {} --species {}'.format(key,
                                                                                                                    "other")
            print("command is:", command)
            print("after printing command")

            os.system(command)
            
            json_results = Redis.blpop("pylmm:results:" + str(temp_uuid), 45*60)
            results = json.loads(json_results[1])
            p_values = [float(result) for result in results['p_values']]
            #print("p_values:", p_values)
            self.dataset.group.markers.add_pvalues(p_values)
            
            self.trait_results[this_trait.name] = []
            for marker in self.dataset.group.markers.markers:
                self.trait_results[this_trait.name].append(marker['lod_score'])

        #print("self.trait_results:", self.trait_results)
            
        chrnames = []
        lodnames = []
        chr_pos = []
        pos = []
        markernames = []
        
        for trait in self.trait_results.keys():
            lodnames.append(trait)
        
        for marker in self.dataset.group.markers.markers:
            if marker['chr'] not in chrnames:
                chrnames.append(marker['chr'])
            chr_pos.append(marker['chr'])
            pos.append(marker['Mb'])
            markernames.append(marker['name'])
            
        self.json_data['chrnames'] = chrnames
        self.json_data['lodnames'] = lodnames
        self.json_data['chr'] = chr_pos
        self.json_data['pos'] = pos
        self.json_data['markernames'] = markernames
        
        for trait in self.trait_results:
            self.json_data[trait] = self.trait_results[trait]
            
        self.js_data = dict(
            json_data = self.json_data
        )
            
        print("self.js_data:", self.js_data)
        
    def identify_empty_samples(self, values):
        no_val_samples = []
        for sample_count, val in enumerate(values):
            if val == "x":
                no_val_samples.append(sample_count)
        return no_val_samples
        
    def trim_genotypes(self, genotype_data, no_value_samples):
        trimmed_genotype_data = []
        for marker in genotype_data:
            new_genotypes = []
            for item_count, genotype in enumerate(marker):
                if item_count in no_value_samples:
                    continue
                try:
                    genotype = float(genotype)
                except ValueError:
                    genotype = np.nan
                    pass
                new_genotypes.append(genotype)
            trimmed_genotype_data.append(new_genotypes)
        return trimmed_genotype_data
            
            