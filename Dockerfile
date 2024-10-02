FROM python:3.8

RUN mkdir /webapp
COPY . /webapp

# installing psql such that it is possible to access the django database
# from within the container
RUN apt-get update && apt-get install -y gnupg2
RUN apt-get install -y postgresql postgresql-contrib

# installing zip to avoid zipping from within python
RUN apt-get install zip

# Lynx
RUN git clone https://github.com/SysMedOs/LipidLynxX.git
RUN cp -r LipidLynxX/lynx/ /usr/local/lib/python3.8/site-packages/

WORKDIR /webapp

RUN pip3 install --no-cache-dir -r requirements.txt

WORKDIR /webapp/lipid_network_project

# improved logging speed
ENV PYTHONBUFFERED 1

#RUN python3 manage.py makemigrations

#RUN python3 manage.py migrate

EXPOSE 7000

ENTRYPOINT ["bash", "startup.sh"]
