import threading
import mysql.connector
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
from ispyb.exception import (ISPyBConnectionException, ISPyBNoResultException,
                             ISPyBRetrieveFailed, ISPyBWriteFailed)
import sshtunnel
import time
import pymysql


class SSHConnector(Connector):
    def __init__(self,
                 user=None,
                 pw=None,
                 host="localhost",
                 db=None,
                 port=3306,
                 reconn_attempts=6,
                 reconn_delay=1,
                 remote=False,
                 ssh_user=None,
                 ssh_password=None,
                 ssh_host=None,
                 conn_inactivity=360,
                 ):
        self.conn_inactivity = conn_inactivity
        self.lock = threading.Lock()
        self.server = None

        if remote:
            creds = {'ssh_host': ssh_host,
                     'ssh_user': ssh_user,
                     'ssh_pass': ssh_password,
                     'db_host': host,
                     'db_port': int(port),
                     'db_user': user,
                     'db_pass': pw,
                     'db_name': db}
            self.remote_connect(**creds)

        else:
            self.connect(user=user, pw=pw, host=host, db=db, port=port, conn_inactivity=conn_inactivity)

    def remote_connect(self, ssh_host, ssh_user, ssh_pass, db_host, db_port, db_user, db_pass, db_name):

        sshtunnel.SSH_TIMEOUT = 5.0
        sshtunnel.TUNNEL_TIMEOUT = 5.0
        self.conn_inactivity = int(self.conn_inactivity)

        self.server = sshtunnel.SSHTunnelForwarder(
            (ssh_host),
            ssh_username=ssh_user, ssh_password=ssh_pass,
            remote_bind_address=(db_host, db_port)
        )

        self.server.daemon_forward_servers = True

        self.server.start()

        self.conn = pymysql.connect(
            user=db_user, password=db_pass,
            host='127.0.0.1', port=self.server.local_bind_port,
            database=db_name,
        )

        if self.conn is not None:
            self.conn.autocommit = True
        else:
            self.server.stop()
            raise ISPyBConnectionException
        self.last_activity_ts = time.time()


