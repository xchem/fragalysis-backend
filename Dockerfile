FROM abradle/rdkit:stable
ENV PYTHONUNBUFFERED 1
RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/
RUN apt-get update -y
RUN pip install -r requirements.txt
# Conver this into a single pip install command
RUN git clone https://github.com/xchem/fragalysis /usr/local/fragalysis
RUN pip install -r /usr/local/fragalysis/requirements.txt
RUN pip install /usr/local/fragalysis
RUN apt-get update -y
RUN apt-get install -y nginx curl
# Copy entrypoint script into the image
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
RUN echo "daemon off;" >> /etc/nginx/nginx.conf
ADD . /code/
RUN mkdir /srv/logs/
