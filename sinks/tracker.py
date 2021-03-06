#!/usr/bin/env python

import fileinput
import sys
import re
import socket
import logging


"""
Simple stream sink to read stdin for incoming metrics of the format
`key|val|timestamp\n`

and tracks the number of unique timers, gauges and counters

Finally sends the gauges to localhost:8125
"""

class Tracker(object):
    GAUGE = "g"

    def __init__(self, host="localhost", port=8125, attempts=3):
        """
        Raises a :class:`ValueError` on bad arguments.

        :Parameters:
            - `host` : The hostname of the statsd server.
            - `port` : The port of the statsd server
            - `attempts` (optional) : The number of re-connect retries before failing.
        """
        # Convert the port to an int since its coming from a configuration file
        port = int(port)
        attempts = int(attempts)

        if port <= 0:
            raise ValueError("Port must be positive!")
        if attempts <= 1:
            raise ValueError("Must have at least 1 attempt!")

        self.addr = (host, port)
        self.attempts = attempts
        self.sock = self._create_socket()
        self.r_counter = re.compile(r'^counters\..*\.sum$')
        self.r_gauge = re.compile(r'^gauges\..*\.sum$')
        self.r_timer = re.compile(r'^timers\..*\.count$')
        self.n_counter = 0
        self.n_gauge = 0
        self.n_timer = 0
        self.n_lines = 0
        self.logger = logging.getLogger("statsite.tracker")

    def measure(self, metric):
        """
        count unique metrics of each type

        :Parameters:
         - `metric` :  statsd metric  "key|value|timestamp" string.
        """
        if not metric:
            return

        # Construct the output
        if metric and metric.count("|") == 2:
            metric = metric.split("|")[0]

        if self.r_gauge.match(metric):
            self.n_gauge += 1
        elif self.r_counter.match(metric):
            self.n_counter += 1
        elif self.r_timer.match(metric):
            self.n_timer += 1

        self.n_lines += 1

    def close(self):
        """
        Closes the connection. The socket will be recreated on the next
        flush.
        """
        self.sock.close()

    def _create_socket(self):
        """Creates a socket and connects to the statsrelay-base"""
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            return sock
        except socket.error:
            self.logger.exception("Failed to connect to statsrelay-base")
            sys.exit(1)

    def flush_metrics(self):
        """Tries to write a string to the socket, reconnecting on any errors"""
        # Write gauges to statsrelay for tracking
        message = list()
        message.append("statsite.timer_s_:{0}|{1}".format(self.n_timer, self.GAUGE))
        message.append("statsite.gauge_s_:{0}|{1}".format(self.n_gauge, self.GAUGE))
        message.append("statsite.counter_s_:{0}|{1}".format(self.n_counter, self.GAUGE))

        self.logger.info("Processing %d metrics" % self.n_lines)

        for _ in range(self.attempts):
            try:
                for item in message:
                    self.sock.sendto(item.encode('utf-8'), self.addr)

                return
            except socket.error:
                self.logger.error('Failed to send')
                self.sock = self._create_socket()

        self.logger.critical("Failed to write metrics to statsrelay. Gave up after %d attempts." % self.attempts)


if __name__ == "__main__":
    # Initialize the logger
    logging.basicConfig()

    # Intialize from our arguments
    tracker = Tracker(*sys.argv[1:])

    # Get all the inputs
    line = 0
    for metric in fileinput.input():
        tracker.measure(metric.rstrip())
        line += 1

    tracker.flush_metrics()

    # Flush
    tracker.close()
