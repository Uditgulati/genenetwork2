import requests
import simplejson as json
from wqflask import app

import os
import string
import resource
import codecs
import requests
import random
import urllib
import datetime as dt

from base import webqtlConfig
from base.webqtlCaseData import webqtlCaseData
from base.data_set import create_dataset
from utility import hmac
from utility.authentication_tools import check_resource_availability
from utility.tools import GN2_BASE_URL
from utility.redis_tools import get_redis_conn, get_resource_id
from utility.tools import GN2_BASE_URL, GN_VERSION
from utility.redis_tools import get_redis_conn
from utility.redis_tools import get_resource_id
from utility.redis_tools import get_resource_info
from utility.benchmark import Bench

from utility.db_tools import escape

from flask import g, request, url_for

from utility.logger import getLogger

logger = getLogger(__name__)

Redis = get_redis_conn()


def create_trait(**kw):
    assert bool(kw.get('dataset')) != bool(
        kw.get('dataset_name')), "Needs dataset ob. or name"

    tb = None
    tc = None
    permitted = True
    dataset = None
    xx1 = 0
    xx2 = 0
    if kw.get('name'):
        if kw.get('dataset_name'):
            if kw.get('dataset_name') != "Temp":
                dataset = create_dataset(kw.get('dataset_name'))
        else:
            dataset = kw.get('dataset')

        if dataset.type != "Temp":
            if dataset.type == 'Publish':
                t1 = dt.datetime.now()
                permissions = check_resource_availability(
                    dataset, kw.get('name'))
                t2 = dt.datetime.now()
                tb = t2 - t1
                xx1 += 1
            else:
                t1 = dt.datetime.now()
                permissions = kw.get('permissions')
                t2 = dt.datetime.now()
                tb = t2 - t1
                xx2 += 1

    if "view" in permissions['data']:
        t1 = dt.datetime.now()
        trait = GeneralTrait(**kw)
        trait_name=kw.get('name')
        get_qtl_info=kw.get('get_qtl_info')

        if dataset.type != "Temp" and get_qtl_info:
            # LRS and its location
            trait.locus = trait.locus_chr = trait.locus_mb = ""
            if dataset.type == 'ProbeSet':
                query = """
                        SELECT 
                            Geno.Chr, Geno.Mb 
                        FROM
                            Geno, Species, ProbeSetXRef, ProbeSet
                        WHERE
                            Species.Name = "{}" AND
                            ProbeSet.Name = "{}" AND
                            ProbeSetXRef.ProbeSetFreezeId = {} AND
                            ProbeSet.Id = ProbeSetXRef.ProbeSetId AND
                            ProbeSetXRef.Locus = Geno.Name AND
                            Geno.SpeciesId = Species.Id
                        """.format(dataset.species, trait_name, dataset.id)
                logger.sql(query)
                result = g.db.execute(query).fetchone()
                if result:
                    trait.locus_chr = result[0]
                    trait.locus_mb = result[1]

            elif dataset.type == 'Publish':
                query = """
                        SELECT
                            Geno.Chr, Geno.Mb
                        FROM
                            Geno, Species, PublishXRef, PublishFreeze
                        WHERE
                            Species.Name = "{}" AND
                            PublishXRef.Id = {} AND
                            PublishFreeze.Id = {} AND
                            Geno.Name = PublishXRef.Locus AND
                            Geno.SpeciesId = Species.Id AND
                            PublishXRef.InbredSetId = PublishFreeze.InbredSetId
                        """ % (dataset.species, trait_name, dataset.id)
                logger.sql(query)
                result = g.db.execute(query).fetchone()
                if result:
                    trait.locus_chr = result[0]
                    trait.locus_mb = result[1]                    

        t2 = dt.datetime.now()
        tc = t2 - t1
        trait.tb = tb
        trait.tc = tc
        trait.xx1 = xx1
        trait.xx2 = xx2

        return trait
    else:
        return None


class GeneralTrait(object):
    """
    Trait class defines a trait in webqtl, can be either Microarray,
    Published phenotype, genotype, or user input trait

    """

    def __init__(self, get_qtl_info=False, get_sample_info=True, **kw):
        # xor assertion
        assert bool(kw.get('dataset')) != bool(
            kw.get('dataset_name')), "Needs dataset ob. or name"
        # Trait ID, ProbeSet ID, Published ID, etc.

    def export_informative(self, include_variance=0):
        """
        export informative sample
        mostly used in qtl regression

        """
        samples = []
        vals = []
        the_vars = []
        sample_aliases = []
        for sample_name, sample_data in list(self.data.items()):
            if sample_data.value is not None:
                if not include_variance or sample_data.variance is not None:
                    samples.append(sample_name)
                    vals.append(sample_data.value)
                    the_vars.append(sample_data.variance)
                    sample_aliases.append(sample_data.name2)
        return samples, vals, the_vars, sample_aliases

    @property
    def description_fmt(self):
        """Return a text formated description"""
        if self.dataset.type == 'ProbeSet':
            if self.description:
                formatted = self.description
                if self.probe_target_description:
                    formatted += "; " + self.probe_target_description
            else:
                formatted = "Not available"
        elif self.dataset.type == 'Publish':
            if self.confidential:
                formatted = self.pre_publication_description
            else:
                formatted = self.post_publication_description
        else:
            formatted = "Not available"
        if isinstance(formatted, bytes):
            formatted = formatted.decode("utf-8")
        return formatted

    @property
    def alias_fmt(self):
        """Return a text formatted alias"""

        alias = 'Not available'
        if getattr(self, "alias", None):
            alias = self.alias.replace(";", " ")
            alias = ", ".join(alias.split())

        return alias

    @property
    def wikidata_alias_fmt(self):
        """Return a text formatted alias"""

        alias = 'Not available'
        if self.symbol:
            human_response = requests.get(
                GN2_BASE_URL + "gn3/gene/aliases/" + self.symbol.upper())
            mouse_response = requests.get(
                GN2_BASE_URL + "gn3/gene/aliases/" + self.symbol.capitalize())
            other_response = requests.get(
                GN2_BASE_URL + "gn3/gene/aliases/" + self.symbol.lower())

            if human_response and mouse_response and other_response:
                alias_list = json.loads(human_response.content) + json.loads(
                    mouse_response.content) + \
                    json.loads(other_response.content)

                filtered_aliases = []
                seen = set()
                for item in alias_list:
                    if item in seen:
                        continue
                    else:
                        filtered_aliases.append(item)
                        seen.add(item)
                alias = "; ".join(filtered_aliases)

        return alias

    @property
    def location_fmt(self):
        """Return a text formatted location

        While we're at it we set self.location in case we need it
        later (do we?)

        """

        if self.chr and self.mb:
            self.location = 'Chr %s @ %s Mb' % (self.chr, self.mb)
        elif self.chr:
            self.location = 'Chr %s @ Unknown position' % (self.chr)
        else:
            self.location = 'Not available'

        fmt = self.location
        # XZ: deal with direction
        if self.strand_probe == '+':
            fmt += (' on the plus strand ')
        elif self.strand_probe == '-':
            fmt += (' on the minus strand ')

        return fmt


def retrieve_sample_data(trait, dataset, samplelist=None):
    if samplelist is None:
        samplelist = []

    if dataset.type == "Temp":
        results = Redis.get(trait.name).split()
    else:
        results = dataset.retrieve_sample_data(trait.name)

    # Todo: is this necessary? If not remove
    trait.data.clear()

    if results:
        if dataset.type == "Temp":
            all_samples_ordered = dataset.group.all_samples_ordered()
            for i, item in enumerate(results):
                try:
                    trait.data[all_samples_ordered[i]] = webqtlCaseData(
                        all_samples_ordered[i], float(item))
                except:
                    pass
        else:
            for item in results:
                name, value, variance, num_cases, name2 = item
                if not samplelist or (samplelist and name in samplelist):
                    # name, value, variance, num_cases)
                    trait.data[name] = webqtlCaseData(*item)
    return trait


@app.route("/trait/get_sample_data")
def get_sample_data():
    params = request.args
    trait = params['trait']
    dataset = params['dataset']

    trait_ob = create_trait(name=trait, dataset_name=dataset)
    if trait_ob:
        trait_dict = {}
        trait_dict['name'] = trait
        trait_dict['db'] = dataset
        trait_dict['type'] = trait_ob.dataset.type
        trait_dict['group'] = trait_ob.dataset.group.name
        trait_dict['tissue'] = trait_ob.dataset.tissue
        trait_dict['species'] = trait_ob.dataset.group.species
        trait_dict['url'] = url_for(
            'show_trait_page', trait_id=trait, dataset=dataset)
        trait_dict['description'] = trait_ob.description_display
        if trait_ob.dataset.type == "ProbeSet":
            trait_dict['symbol'] = trait_ob.symbol
            trait_dict['location'] = trait_ob.location_repr
        elif trait_ob.dataset.type == "Publish":
            if trait_ob.pubmed_id:
                trait_dict['pubmed_link'] = trait_ob.pubmed_link
            trait_dict['pubmed_text'] = trait_ob.pubmed_text

        return json.dumps([trait_dict, {key: value.value for
                                        key, value in list(
                                            trait_ob.data.items())}])
    else:
        return None


def jsonable(trait):
    """Return a dict suitable for using as json

    Actual turning into json doesn't happen here though"""

    dataset = create_dataset(dataset_name=trait.dataset.name,
                             dataset_type=trait.dataset.type,
                             group_name=trait.dataset.group.name)

    if dataset.type == "ProbeSet":
        return dict(name=trait.name,
                    symbol=trait.symbol,
                    dataset=dataset.name,
                    dataset_name=dataset.shortname,
                    description=trait.description_display,
                    mean=trait.mean,
                    location=trait.location_repr,
                    lrs_score=trait.LRS_score_repr,
                    lrs_location=trait.LRS_location_repr,
                    additive=trait.additive
                    )
    elif dataset.type == "Publish":
        if trait.pubmed_id:
            return dict(name=trait.name,
                        dataset=dataset.name,
                        dataset_name=dataset.shortname,
                        description=trait.description_display,
                        abbreviation=trait.abbreviation,
                        authors=trait.authors,
                        pubmed_text=trait.pubmed_text,
                        pubmed_link=trait.pubmed_link,
                        lrs_score=trait.LRS_score_repr,
                        lrs_location=trait.LRS_location_repr,
                        additive=trait.additive
                        )
        else:
            return dict(name=trait.name,
                        dataset=dataset.name,
                        dataset_name=dataset.shortname,
                        description=trait.description_display,
                        abbreviation=trait.abbreviation,
                        authors=trait.authors,
                        pubmed_text=trait.pubmed_text,
                        lrs_score=trait.LRS_score_repr,
                        lrs_location=trait.LRS_location_repr,
                        additive=trait.additive
                        )
    elif dataset.type == "Geno":
        return dict(name=trait.name,
                    dataset=dataset.name,
                    dataset_name=dataset.shortname,
                    location=trait.location_repr
                    )
    else:
        return dict()


def jsonable_table_row(trait, dataset_name, index):
    """Return a list suitable for json and intended to be displayed in a table

    Actual turning into json doesn't happen here though"""

    dataset = create_dataset(dataset_name)

    if dataset.type == "ProbeSet":
        if trait.mean == "":
            mean = "N/A"
        else:
            mean = "%.3f" % round(float(trait.mean), 2)
        if trait.additive == "":
            additive = "N/A"
        else:
            additive = "%.3f" % round(float(trait.additive), 2)
        return ['<input type="checkbox" name="searchResult" class="checkbox trait_checkbox" value="' + hmac.data_hmac('{}:{}'.format(str(trait.name), dataset.name)) + '">',
                index,
                '<a href="/show_trait?trait_id=' +
                str(trait.name)+'&dataset='+dataset.name +
                '">'+str(trait.name)+'</a>',
                trait.symbol,
                trait.description_display,
                trait.location_repr,
                mean,
                trait.LRS_score_repr,
                trait.LRS_location_repr,
                additive]
    elif dataset.type == "Publish":
        if trait.additive == "":
            additive = "N/A"
        else:
            additive = "%.2f" % round(float(trait.additive), 2)
        if trait.pubmed_id:
            return ['<input type="checkbox" name="searchResult" class="checkbox trait_checkbox" value="' + hmac.data_hmac('{}:{}'.format(str(trait.name), dataset.name)) + '">',
                    index,
                    '<a href="/show_trait?trait_id=' +
                    str(trait.name)+'&dataset='+dataset.name +
                    '">'+str(trait.name)+'</a>',
                    trait.description_display,
                    trait.authors,
                    '<a href="' + trait.pubmed_link + '">' + trait.pubmed_text + '</href>',
                    trait.LRS_score_repr,
                    trait.LRS_location_repr,
                    additive]
        else:
            return ['<input type="checkbox" name="searchResult" class="checkbox trait_checkbox" value="' + hmac.data_hmac('{}:{}'.format(str(trait.name), dataset.name)) + '">',
                    index,
                    '<a href="/show_trait?trait_id=' +
                    str(trait.name)+'&dataset='+dataset.name +
                    '">'+str(trait.name)+'</a>',
                    trait.description_display,
                    trait.authors,
                    trait.pubmed_text,
                    trait.LRS_score_repr,
                    trait.LRS_location_repr,
                    additive]
    elif dataset.type == "Geno":
        return ['<input type="checkbox" name="searchResult" class="checkbox trait_checkbox" value="' + hmac.data_hmac('{}:{}'.format(str(trait.name), dataset.name)) + '">',
                index,
                '<a href="/show_trait?trait_id=' +
                str(trait.name)+'&dataset='+dataset.name +
                '">'+str(trait.name)+'</a>',
                trait.location_repr]
    else:
        return dict()


def retrieve_trait_info(trait, dataset, trait_name=None, trait_cellid=None, get_qtl_info=False):
    assert dataset, "Dataset doesn't exist"

    if get_qtl_info:
        # LRS and its location
        trait.locus = trait.locus_chr = trait.locus_mb = ""
        if dataset.type == 'ProbeSet' and not trait_cellid:
            trait.mean = ""
            query = """
                    SELECT
                            ProbeSetXRef.Locus
                    FROM
                            ProbeSetXRef, ProbeSet
                    WHERE
                            ProbeSetXRef.ProbeSetId = ProbeSet.Id AND
                            ProbeSet.Name = "{}" AND
                            ProbeSetXRef.ProbeSetFreezeId ={}
                    """.format(trait_name, dataset.id)
            logger.sql(query)
            trait_qtl = g.db.execute(query).fetchone()

        if dataset.type == 'Publish':
            query = """
                    SELECT
                            PublishXRef.Locus
                    FROM
                            PublishXRef, PublishFreeze
                    WHERE
                            PublishXRef.Id = %s AND
                            PublishXRef.InbredSetId = PublishFreeze.InbredSetId AND
                            PublishFreeze.Id =%s
                    """ % (trait_name, dataset.id)
            logger.sql(query)
            trait_qtl = g.db.execute(query).fetchone()

        if trait_qtl:
            trait.locus = trait_qtl[0]
            if trait.locus:
                query = """
                    select Geno.Chr, Geno.Mb from Geno, Species
                    where Species.Name = '{}' and
                    Geno.Name = '{}' and
                    Geno.SpeciesId = Species.Id
                    """.format(dataset.group.species, trait.locus)
                logger.sql(query)
                result = g.db.execute(query).fetchone()
                if result:
                    trait.locus_chr = result[0]
                    trait.locus_mb = result[1]
    return trait
