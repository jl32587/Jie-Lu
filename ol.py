#!/usr/bin/python
def overlap(a,b,percentage = False):
    '''Determines whether two lists of coordinates overlap. a and b are lists of genome coordinates[[start1,end1],[start2,end2]...[startn,endn]] that have already been sorted by start and end.
    The function returns a list of lists that has the same index as b.'''
    # for example:
    # a = [[19,38],[21,50],[200,300],[250,400]]
    # b = [[18,39],[50,80],[200,300]]
    # overlap(a,b) returns a list of lists ol:
    # [[0, 1],[],[2, 3]]
    total_len = sum([i[1]-i[0]+1 for i in a])
    
    m = len(a)
    n = len(b)
    a.append([(),()])
    b.append([(),()])                   # add infinity to the end of each                                                 array so that the last element of                                               either a and b has something to                                                 compare to.
    i = 0
    j = 0    
    ol = [[] for u in range(n+1)]
    ol_len = 0
   
    
    while i <= m and j <= n:      
                                        # if a[i] and b[j] do not overlap
        if a[i][1] < b[j][0]:           # if a[i] is on the left of b[j]
            i += 1                      # move to the next element of a
        elif a[i][0] > b[j][1]:         # if a[i] is on the right of b[j]
            j += 1                      # move to the next element of b
            if b[j][0] < b[j-1][1]:     # trace back all a if current b overlaps                                          with the previous b
               i =  traceback(b[j][0],i,a)
           
                                        # if a[i] and b[j] overlap
        else:
            ol[j].append(i)             # append the index of the element in                                              list a which overlaps with b[j] to                                              the result list ol[j], which means                                              that b[j] overlaps with a[i].
            if percentage == True:
                ol_len += overlap_len(a[i],b[j])
                
            if a[i][1] > b[j][1]:       # if the end of a[i] extend beyond the                                            end of b[j], it could also overlap                                              with the following elements of b.
                for y in range(j+1, n): # loop through all the other elements                                             in b to check whether a[i] overlaps                                             with the rest elements of b.
                    if a[i][1] > b[y][0]:
                        ol[y].append(i)
                        
                        if percentage == True:
                            ol_len += overlap_len(a[i],b[y])
                            
                    else:               # until a[i] does not overlap with any                                            other elements of b.
                        break           
            i += 1                      # move to the next element of a
            
    del ol[n]
    del a[m]
    del b[n]
    if percentage == False:
        return ol
    else:
        return float(ol_len)/total_len

def traceback(b,i,a):
    while i > 0:
        if a[i][1] > b:
            i -= 1
        else:
            break
    return i

def overlap_dot(a,b):
    '''Determines whether two lists of coordinates overlap. a and b are lists of genome coordinates[[start1,end1],[start2,end2]...[startn,endn]] that have already been sorted by start and end.
    The function returns a list of lists that has the same index as b.'''

    # for example:
    # a = [19,21,50,200,300,400]
    # b = [[18,39],[50,80],[200,300]]
    # overlap(a,b) returns a list of lists ol:
 
    m = len(a)
    n = len(b)
    a.append([])
    b.append([(),()])                        # add infinity to the end of each                                                 array so that the last element of                                               either a and b has something to                                                 compare to.
    i = 0
    j = 0    
    ol = [[] for u in range(n+1)]
    
    while i <= m and j <= n:      
                                        # if a[i] and b[j] do not overlap
        if a[i] < b[j][0]:              # if a[i] is on the left of b[j]
            i += 1                      # move to the next element of a
        elif a[i] > b[j][1]:            # if a[i] is on the right of b[j]
            j += 1                      # move to the next element of b
            if b[j][0] < b[j-1][1]:     
               i =  traceback_dot(b[j][0],i,a)
           
                                        # if a[i] and b[j] overlap
        else:
            ol[j].append(i)             # append the index of the element in                                              list a which overlaps with b[j] to                                              the result list ol[j], which means                                              that b[j] overlaps with a[i].       
            i += 1                      # move to the next element of a
            
    del ol[n]
    del a[m]
    del b[n]
    return ol

def traceback_dot(b,i,a):
    while i > 0:
        if a[i] > b:
            i -= 1
        else:
            break
    return i

def overlap_len(a,b):
    try:
        if a[0] <= b[0]:
            if a[1] <= b[1]:
                return a[1] - b[0] + 1
            else:
                return b[1] - b[0] + 1
        else:
            if a[1] <= b[1]:
                return a[1] - a[0] + 1
            else:
                return b[1] - a[0] + 1
    except TypeError:
        return 0


## a = [[19,38],[21,50],[200,300],[250,400]]
## b = [[18,39],[50,80],[200,300]]
## print overlap(a,b,percentage=True)
