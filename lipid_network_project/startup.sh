python3 manage.py makemigrations lipid_network
python3 manage.py makemigrations
python3 manage.py migrate lipid_network
python3 manage.py migrate
python3 -u /webapp/lipid_network_project/manage.py process_tasks &
#python3 /root/webapp/lipid_network_project/manage.py runserver 0.0.0.0:7000
gunicorn lipid_network_project.wsgi:application --bind 0.0.0.0:7000 --timeout 600
