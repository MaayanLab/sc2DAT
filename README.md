# MULTIOMICS2TARGETS: 
### Expression2Kinases & Target Ranger

Currently hosted at: https://multiomics2targets.maayanlab.cloud/


## Getting Started
To run in development:
```bash
# prepare .env file & review
cp .env.example .env
# create and activate python3.8 virtual environment
python3.8 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
# install relevant R dependencies
R -e "source('setup.R')"
# run in development
appyter main.ipynb
```



