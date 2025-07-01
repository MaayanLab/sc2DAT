from flask import Blueprint, make_response
import os
import s3fs
from appyter.ext.fsspec.core import url_to_chroot_fs
from dotenv import load_dotenv
load_dotenv()
blueprint = Blueprint('blueprint', __name__)

@blueprint.route('numruns')
def numruns():
  fs = url_to_chroot_fs(os.environ.get('APPYTER_DATA_DIR'))
  n_runs = len(fs.ls('output'))
  response = make_response({'num_runs': n_runs})
  response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
  response.headers['Pragma'] = 'no-cache'
  response.headers['Expires'] = '0'
  return response