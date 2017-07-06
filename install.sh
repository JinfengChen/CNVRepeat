virtualenv env
source env/bin/activate
export PYTHONPATH=`pwd`/env/lib/python2.7/site-packages
pip install --upgrade pip
pip install .

#./env/bin/CNVRepeat
#./env/bin/CNVRepeat --config config.json
#deactivate env
