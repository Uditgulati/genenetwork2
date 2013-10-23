from __future__ import print_function, division, absolute_import


import os
import hashlib
import datetime
import time

import uuid
import hashlib
import hmac
import base64

import urlparse

import simplejson as json

from sqlalchemy import orm

#from redis import StrictRedis
import redis
Redis = redis.StrictRedis()


from flask import (Flask, g, render_template, url_for, request, make_response,
                   redirect, flash)

from wqflask import app


from pprint import pformat as pf


from wqflask.database import db_session

from wqflask import model

from utility import Bunch, Struct

from wqflask import user_manager






@app.route("/collections/add")
def collections_add():
    return render_template("collections/add.html", traits=request.args['traits'])


@app.route("/collections/new")
def collections_new():
    uc = model.UserCollection()
    uc.name = request.args['new_collection']
    print("user_session:", g.user_session.__dict__)
    uc.user = g.user_session.record['user_id']
    unprocessed_traits = request.args['traits']
    print("unprocessed_traits are:", unprocessed_traits)
    unprocessed_traits = unprocessed_traits.split(",")
    traits = set()
    for trait in unprocessed_traits:
        data, _separator, hmac = trait.rpartition(':')
        data = data.strip()
        print("data is:", data)
        print("hmac is:", hmac)
        assert hmac==user_manager.actual_hmac_creation(data), "Data tampering?"
        traits.add(str(data))

    uc.members = json.dumps(list(traits))
    print("traits are:", traits)

    db_session.add(uc)
    db_session.commit()



    return "Created: " + uc.name
