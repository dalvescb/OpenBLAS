#!/usr/bin/env python
file = open('objlist_tmp', 'r')
Lines = file.readlines()

path=""
for line in Lines:
  if line[0]=='/':
    path=line
  else:
    for f in line.split(' '):
      print(path.strip() + "/" + f.strip())
