FROM xchem/rdkit-python-debian:Release_2017_09_2
ENV PYTHONUNBUFFERED 1
ADD . /code/
WORKDIR /code
USER root
RUN apt-get update -y
RUN apt-get install -y nginx curl git
RUN git clone https://github.com/xchem/fragalysis /usr/local/fragalysis
RUN pip install -r requirements.txt
# Conver this into a single pip install command
RUN pip install -r /usr/local/fragalysis/requirements.txt
RUN pip install /usr/local/fragalysis
# Copy entrypoint script into the image
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
RUN echo "daemon off;" >> /etc/nginx/nginx.conf
RUN mkdir /srv/logs/