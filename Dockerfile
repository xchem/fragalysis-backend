FROM informaticsmatters/rdkit
ENV PYTHONUNBUFFERED 1
RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/
RUN pip install -r requirements.txt
RUN git clone https://github.com/xchem/fragalysis /usr/local/fragalysis
RUN pip install -r /usr/local/fragalysis/requirements.txt
RUN pip install /usr/local/fragalysis
RUN apt-get update -y
RUN apt-get install -y nginx curl
RUN curl -sL https://deb.nodesource.com/setup_6.x | bash -
RUN apt-get install -y nodejs
# Copy entrypoint script into the image
COPY django_nginx.conf /etc/nginx/sites-available/default.conf
RUN ln -s /etc/nginx/sites-available/default.conf /etc/nginx/sites-enabled
RUN echo "daemon off;" >> /etc/nginx/nginx.conf
ADD . /code/
RUN mkdir /srv/logs/
