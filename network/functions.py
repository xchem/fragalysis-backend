import psycopg2
def get_conn():
    conn = psycopg2.connect(database='dsi',port=5433,host='52.91.71.182',user="postgres")
    return conn