#!/bin/bash
gcc dlinsol.c -I/usr/include/superlu/ -L/usr/lib/x86_64-linux-gnu/ -lsuperlu -g
