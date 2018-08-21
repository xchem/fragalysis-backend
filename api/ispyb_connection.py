from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
import os


def get_conn():
    credentials = {
        "user": os.environ["ISPYB_USER"],
        "pw": os.environ["ISPYB_PASSWORD"],
        "host": os.environ["ISPYB_HOST"],
        "port": os.environ["ISPYB_PORT"],
        "db": os.environ["ISPYB_DB"],
        "conn_inactivity": 360,
    }
    conn = Connector(**credentials)
    return conn
