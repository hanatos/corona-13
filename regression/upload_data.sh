#!/bin/bash

tar cvf data.tar.gz $(cat datafiles.txt) && scp data.tar.gz dreggn:/srv/web/jo.dreggn.org/corona-regression-data.tar.gz && rm data.tar.gz
