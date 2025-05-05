# sc2DAT: 
### Single Cells to Drugs and Targets


## Getting Started
To run in development:
```bash
# prepare .env file & review
cp .env.example .env
# create and activate python3.8 (or python3.9) virtual environment
python3.9 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
# install relevant R dependencies
R -e "source('setup.R')"
# run in development
appyter main.ipynb
```

## Running locally
To run the SC2Targets app locally first ensure that [Docker](https://www.docker.com/) is installed, and then run the following command:
```bash
docker run --device /dev/fuse --cap-add SYS_ADMIN --security-opt apparmor:unconfined -p 5000:5000 -it maayanlab/sc2dat:0.0.18
```

## Preparing Additional Single Cell References
To prepare additional scRNA-seq reference backgrounds:
- `prepareRef.py` - If downloaded in an h5ad format utilize this script
- `prepareRef.r` - If downloaded in an RDS format utilize this script

Once you have prepared an appropriate .rds reference file, it can be added to the reference list in main.ipynb and can be selected from the input interface after running the app locally.
