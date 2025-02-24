#!/bin/bash

while true; do
    free -m | awk '/^Mem:/ {print strftime("%Y-%m-%d %H:%M:%S"), "Total Used:", $3, "MB"}' >> ram_usage.log
    sleep 10
done