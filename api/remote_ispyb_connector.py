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
PYMYSQL_CONNECT_TIMEOUT_S = 3
PYMYSQL_READ_TIMEOUT_S = 3
PYMYSQL_WRITE_TIMEOUT_S = 10
# MySQL DB connection attempts.
# An attempt to cope with intermittent OperationalError exceptions
# that are seen to occur at "busy times". See m2ms-1403.
PYMYSQL_OE_RECONNECT_ATTEMPTS = 5
PYMYSQL_EXCEPTION_RECONNECT_DELAY_S = 1


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
            logger.debug("Creating remote connector: %s", creds)
            self.remote_connect(**creds)
            logger.debug(
                "Started remote ssh_host=%s ssh_user=%s local_bind_port=%s",
                ssh_host,
                ssh_user,
                self.server.local_bind_port,
            )
        else:
            logger.debug("Creating connector")
            self.connect(
                user=user,
                pw=pw,
                host=host,
                db=db,
                port=port,
                conn_inactivity=conn_inactivity,
            )
            logger.debug("Started direct host=%s user=%s port=%s", host, user, port)

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
        sshtunnel.DEFAULT_LOGLEVEL = logging.ERROR
        self.conn_inactivity = int(self.conn_inactivity)

        if ssh_pkey:
            logger.debug(
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
            logger.debug(
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
        logger.debug('Created SSHTunnelForwarder')

        # stops hanging connections in transport
        self.server.daemon_forward_servers = True
        self.server.daemon_transport = True

        logger.debug('Starting SSH server...')
        self.server.start()
        logger.debug('Started SSH server')

        # Try to connect to the database
        # a number of times (because it is known to fail)
        # before giving up...
        connect_attempts = 0
        self.conn = None
        while self.conn is None and connect_attempts < PYMYSQL_OE_RECONNECT_ATTEMPTS:
            try:
                self.conn = pymysql.connect(
                    user=db_user,
                    password=db_pass,
                    host='127.0.0.1',
                    port=self.server.local_bind_port,
                    database=db_name,
                    connect_timeout=PYMYSQL_CONNECT_TIMEOUT_S,
                    read_timeout=PYMYSQL_READ_TIMEOUT_S,
                    write_timeout=PYMYSQL_WRITE_TIMEOUT_S,
                )
            except OperationalError as oe_e:
                if connect_attempts == 0:
                    # So we only log our connection attempts once
                    # an error has occurred - to avoid flooding the log
                    logger.info(
                        'Connecting to MySQL database (db_user=%s db_name=%s)...',
                        db_user,
                        db_name,
                    )
                logger.warning('%s', repr(oe_e))
                connect_attempts += 1
                time.sleep(PYMYSQL_EXCEPTION_RECONNECT_DELAY_S)
            except Exception as e:
                if connect_attempts == 0:
                    # So we only log our connection attempts once
                    # an error has occurred - to avoid flooding the log
                    logger.info(
                        'Connecting to MySQL database (db_user=%s db_name=%s)...',
                        db_user,
                        db_name,
                    )
                logger.warning('Unexpected %s', repr(e))
                connect_attempts += 1
                time.sleep(PYMYSQL_EXCEPTION_RECONNECT_DELAY_S)

        if self.conn is not None:
            if connect_attempts > 0:
                logger.info('Connected')
            self.conn.autocommit = True
        else:
            if connect_attempts > 0:
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