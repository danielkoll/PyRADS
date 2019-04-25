from __future__ import division, print_function
from copy import copy               # to copy lists of objects...

# ==================================================================================
# Retrieve object(s) from obj list that match obj.attribute == value
# or function(obj.attribute,value).
# Either return object or list.
# Usages:
# sublist = get_objects_from_list(masterlist,"name","Dry_50")
# sublist = get_objects_from_list(masterlist,"omega",5,comparefn=largerthan)
# (note: comparefn(x,y) must return a boolean!)

def get_objects_from_list(obj_list,attribute,value,comparefn=None):
    if comparefn is None:
        x = [obj for obj in obj_list if getattr(obj,attribute)==value]
    else:
        x = [obj for obj in obj_list if comparefn(getattr(obj,attribute),value)]

    if len(x)==1:
        x = x[0]
    return x
