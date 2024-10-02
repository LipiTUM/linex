# LINEX - The Lipid Network EXplorer

Flexible Django webapp, to compute and visualize lipid metabolic networks with statistical properties.


**Publications**: 

* LINEX version 2: [Rose and Köhler et al. **Lipid network and moiety analysis for revealing enzymatic dysregulation and mechanistic alterations from lipidomics data**. 
 _Briefings in Bioinformatics_, 2023, bbac572](https://doi.org/10.1093/bib/bbac572)

* LINEX: [Köhler and Rose et al. **Investigating Global Lipidome Alterations with the Lipid Network Explorer.**
            _Metabolites_ 2021, 11, 488](https://doi.org/10.3390/metabo11080488)

We recommend running LINEX in a docker environment, to avoid compatibility issues.
Below are instructions to run in with and without Docker:

## Deployment with Docker

[docker](https://docs.docker.com/engine/install/) and [docker-compose](https://docs.docker.com/compose/install/) need to be installed.

Depending on your installation you might need to execute the docker commands with sudo (or as administrator)

To start LINEX, run the following command:
```
docker-compose up -d
```
Open your browser at: [http://127.0.0.1:8084](http://127.0.0.1:8084)

To shut down the server run:
```
docker-compose down
```
(Changing the code of the django app requires rebuilding of the web container.)

## Deployment without Docker

### 1. Requirements
* python: >= 3.8
* LipidLynxX (version from October 2020)
    * download lynx folder from [github](https://github.com/SysMedOs/LipidLynxX/)
    * move folder into the site-packages directory of your current environment
        * if you need to find the site-packages dir run (assuming numpy is installed):
            ```shell
            python3 -c "import numpy as np; import os; print(os.path.dirname(os.path.dirname(np.__file__)))"
            ```
### 2. Setup

#### 2.1. Installation
```shell
cd lipid_network_project  # directory containing manage.py
python3 manage.py makemigrations lipid_network
python3 manage.py makemigrations
python3 manage.py migrate lipid_network
python3 manage.py migrate
```

#### 2.2. Running
__Important:__ runserver and process_tasks need to be run in different threads!
```shell
python3 manage.py runserver 7000 &
python3 manage.py process_tasks
```

#### 2.3. Admin Account
To view the content of the database and saved userdata, you need an admin account, which you can create with:
```shell
python3 manage.py createsuperuser
```

### 3. View application
Open your browser at: [http://127.0.0.1:7000](http://127.0.0.1:7000)

## License

LINEX is free software published under the AGPLv3 license.

![AGPLv3 logo](https://www.gnu.org/graphics/agplv3-with-text-162x68.png)
