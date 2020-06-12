#!/usr/bin/env bash
docker rm stochsscontainer_ecoli_signaling
#IMAGE_NAME='briandrawert/stochss-launcher:1.9'
IMAGE_NAME='briandrawert/stochss_bacterial_signaling'

bash -c "sleep 6;xdg-open http://localhost:9999" &

docker run -it -p 9999:9999 --name=stochsscontainer_ecoli_signaling --volume "`pwd`":/working $IMAGE_NAME bash -c "cd /working; C_INCLUDE_PATH=/usr/lib/openmpi/include PYTHONPATH=/pyurdme jupyter notebook --NotebookApp.open_browser=False --ip=0.0.0.0 --port=9999 --config=/stochss-master/jupyter_profile/jupyter.config"
