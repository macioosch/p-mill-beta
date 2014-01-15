#!/usr/bin/env bash
echo $(grep '\\' $@ | tr '\n' ' ' | sed 's/\\//g;s/  */ /g')
