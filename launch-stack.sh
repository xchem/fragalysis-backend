#!/bin/bash
/bin/bash /code/makemigrations.sh
echo "Creating superuser..."
# Automatically create the superuser...
script="
from django.contrib.auth.models import User;
username = '$WEB_DJANGO_SUPERUSER_NAME';
password = '$WEB_DJANGO_SUPERUSER_PASSWORD';
email = '$WEB_DJANGO_SUPERUSER_EMAIL';
if User.objects.filter(username=username).count()==0:
    User.objects.create_superuser(username, email, password);
    print('Superuser created.');
else:
    print('Superuser creation skipped.');
"
printf "$script" | python manage.py shell

# Prepare log files and start outputting logs to stdout
touch /srv/logs/gunicorn.log
touch /srv/logs/access.log
tail -n 0 -f /srv/logs/*.log &
# Start the NPM build
echo "Starting..."
# TODO set up flag that cun this in dev / prod mode
cd /code/frontend && npm run dev &
# Start Gunicorn processes
echo "Starting Gunicorn...."
gunicorn fragalysis.wsgi:application \
    --daemon \
    --name fragalysis \
    --bind unix:django_app.sock \
    --timeout 60 \
    --workers 3 \
    --log-level=debug \
    --log-file=/srv/logs/gunicorn.log \
    --access-logfile=/srv/logs/access.log
echo "Testing nginx config..."
nginx -t
echo "Running nginx..."
nginx