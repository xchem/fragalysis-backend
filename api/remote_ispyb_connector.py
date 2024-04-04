import logging
import threading
import time
import traceback

import pymysql
import sshtunnel
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
from ispyb.exception import (
    ISPyBConnectionException,
    ISPyBNoResultException,
    ISPyBRetrieveFailed,
)
from pymysql.err import OperationalError

logger: logging.Logger = logging.getLogger(__name__)

# Timeout to allow the pymysql.connect() method to connect to the DB.
# The default, if not specified, is 10 seconds.
PYMYSQL_CONNECT_TIMEOUT_S = 5
# MySQL DB connection attempts.
# An attempt to cope with intermittent OperationalError exceptions
# that are seen to occur at "busy times". See m2ms-1403.
PYMYSQL_OE_RECONNECT_ATTEMPTS = 3
PYMYSQL_OE_RECONNECT_DELAY_S = 1


class SSHConnector(Connector):
    def __init__(
        self,
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
        ssh_private_key_filename=None,
        ssh_host=None,
        conn_inactivity=360,
    ):
        # Unused arguments
        del reconn_attempts, reconn_delay

        self.conn_inactivity = conn_inactivity
        self.lock = threading.Lock()
        self.conn = None
        self.server = None
        self.last_activity_ts = None

        if remote:
            creds = {
                'ssh_host': ssh_host,
                'ssh_user': ssh_user,
                'ssh_pass': ssh_password,
                'ssh_pkey': ssh_private_key_filename,
                'db_host': host,
                'db_port': int(port),
                'db_user': user,
                'db_pass': pw,
                'db_name': db,
            }
            logger.info("Creating remote connector: %s", creds)
            self.remote_connect(**creds)
            logger.info(
                "Started remote ssh_host=%s ssh_user=%s local_bind_port=%s",
                ssh_host,
                ssh_user,
                self.server.local_bind_port,
            )
        else:
            logger.info("Creating connector")
            self.connect(
                user=user,
                pw=pw,
                host=host,
                db=db,
                port=port,
                conn_inactivity=conn_inactivity,
            )
            logger.info("Started direct host=%s user=%s port=%s", host, user, port)

    def remote_connect(
        self,
        ssh_host,
        ssh_user,
        ssh_pass,
        ssh_pkey,
        db_host,
        db_port,
        db_user,
        db_pass,
        db_name,
    ):
        sshtunnel.SSH_TIMEOUT = 5.0
        sshtunnel.TUNNEL_TIMEOUT = 5.0
        sshtunnel.DEFAULT_LOGLEVEL = logging.DEBUG
        self.conn_inactivity = int(self.conn_inactivity)

        if ssh_pkey:
            logger.info(
                'Creating SSHTunnelForwarder (with SSH Key) host=%s user=%s',
                ssh_host,
                ssh_user,
            )
            self.server = sshtunnel.SSHTunnelForwarder(
                (ssh_host),
                ssh_username=ssh_user,
                ssh_pkey=ssh_pkey,
                remote_bind_address=(db_host, db_port),
            )
        else:
            logger.info(
                'Creating SSHTunnelForwarder (with password) host=%s user=%s',
                ssh_host,
                ssh_user,
            )
            self.server = sshtunnel.SSHTunnelForwarder(
                (ssh_host),
                ssh_username=ssh_user,
                ssh_password=ssh_pass,
                remote_bind_address=(db_host, db_port),
            )
        logger.info('Created SSHTunnelForwarder')

        # stops hanging connections in transport
        self.server.daemon_forward_servers = True
        self.server.daemon_transport = True

        logger.info('Starting SSH server...')
        self.server.start()
        logger.info('Started SSH server')

        logger.info('Connecting to ISPyB (db_user=%s db_name=%s)...', db_user, db_name)
        # Try to connect to the database
        # a number of times (because it is known to fail)
        # before giving up...
        connect_attempts = 0
        stop_connecting = False
        self.conn = None
        while not stop_connecting and connect_attempts < PYMYSQL_OE_RECONNECT_ATTEMPTS:
            try:
                self.conn = pymysql.connect(
                    user=db_user,
                    password=db_pass,
                    host='127.0.0.1',
                    port=self.server.local_bind_port,
                    database=db_name,
                    connect_timeout=PYMYSQL_CONNECT_TIMEOUT_S,
                )
            except OperationalError as oe_e:
                logger.info(
                    'OperationalError(%d): %s',
                    oe_e.args[0],
                    repr(oe_e),
                )
                connect_attempts += 1
                time.sleep(PYMYSQL_OE_RECONNECT_DELAY_S)
            except Exception as e:
                logger.info(
                    'Unexpected Exception(%d): %s',
                    e.args[0],
                    repr(e),
                )
                stop_connecting = True

        if self.conn is not None:
            logger.info('Connected')
            self.conn.autocommit = True
        else:
            logger.info('Failed to connect')
            self.server.stop()
            raise ISPyBConnectionException
        self.last_activity_ts = time.time()

    def create_cursor(self):
        if time.time() - self.last_activity_ts > self.conn_inactivity:
            # re-connect:
            self.connect(self.user, self.pw, self.host, self.db, self.port)
        self.last_activity_ts = time.time()
        if self.conn is None:
            raise ISPyBConnectionException

        cursor = self.conn.cursor(pymysql.cursors.DictCursor)
        if cursor is None:
            raise ISPyBConnectionException
        return cursor

    def call_sp_retrieve(self, procname, args):
        with self.lock:
            cursor = self.create_cursor()
            try:
                cursor.callproc(procname=procname, args=args)
            except self.conn.DataError as e:
                raise ISPyBRetrieveFailed(
                    "DataError({0}): {1}".format(e.errno, traceback.format_exc())
                ) from e

            result = cursor.fetchall()

            cursor.close()
        if result == []:
            raise ISPyBNoResultException
        return result

    def stop(self):
        if self.server is not None:
            self.server.stop()
        self.server = None
        self.conn = None
        self.last_activity_ts = None
        logger.debug("Server stopped")
