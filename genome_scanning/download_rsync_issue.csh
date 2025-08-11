# 15 retries
set dfile1=$1

set dfile2=$2

set domain=$3

set assembly=$4

set filet1=$5

set filet2=$6


rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/

set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`

if ($dcheck == "0") then
 sleep 60s
 rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
 set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
 if ($dcheck == "0") then
  sleep 60s
  rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
  set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
  if ($dcheck == "0") then
   sleep 60s
   rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
   set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
   if ($dcheck == "0") then
    sleep 60s
    rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
    set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
    if ($dcheck == "0") then
     sleep 60s
     rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
     set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
     if ($dcheck == "0") then
      sleep 60s
      rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
      set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
      if ($dcheck == "0") then
       sleep 60s
       rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
       set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
       if ($dcheck == "0") then
        sleep 60s
        rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
        set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
        if ($dcheck == "0") then
         sleep 60s
         rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
         set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
         if ($dcheck == "0") then
          sleep 60s
          rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
          set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
          if ($dcheck == "0") then
           sleep 60s
           rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
           set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
           if ($dcheck == "0") then
            sleep 60s
            rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
            set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
            if ($dcheck == "0") then
             sleep 60s
             rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
             set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
             if ($dcheck == "0") then
              sleep 60s
              rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
              set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
              if ($dcheck == "0") then
               sleep 60s
               rsync --copy-links --times --verbose $dfile1  $domain/"$assembly"/
               set dcheck=`ls $domain/"$assembly" | grep $filet1 | wc -l`
              endif
             endif
            endif
           endif
          endif
         endif
        endif
       endif
      endif
     endif
    endif
   endif
  endif
 endif
endif



rsync --copy-links --times --verbose $dfile2 $domain/"$assembly"/

set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`

if ($dcheck2 == "0") then
 sleep 60s
 rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
 set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
 if ($dcheck2 == "0") then
  sleep 60s
  rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
  set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
  if ($dcheck2 == "0") then
   sleep 60s
   rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
   set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
   if ($dcheck2 == "0") then
    sleep 60s
    rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
    set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
    if ($dcheck2 == "0") then
     sleep 60s
     rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
     set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
     if ($dcheck2 == "0") then
      sleep 60s
      rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
      set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
      if ($dcheck2 == "0") then
       sleep 60s
       rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
       set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
       if ($dcheck2 == "0") then
        sleep 60s
        rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
        set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
        if ($dcheck2 == "0") then
         sleep 60s
         rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
         set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
         if ($dcheck2 == "0") then
          sleep 60s
          rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
          set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
          if ($dcheck2 == "0") then
           sleep 60s
           rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
           set dcheck=`ls $domain/"$assembly" | grep $filet2 | wc -l`
           if ($dcheck2 == "0") then
            sleep 60s
            rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
            set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
            if ($dcheck2 == "0") then
             sleep 60s
             rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
             set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
             if ($dcheck2 == "0") then
              sleep 60s
              rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
              set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
              if ($dcheck2 == "0") then
               sleep 60s
               rsync --copy-links --times --verbose $dfile2  $domain/"$assembly"/
               set dcheck2=`ls $domain/"$assembly" | grep $filet2 | wc -l`
              endif
             endif
            endif
           endif
          endif
         endif
        endif
       endif
      endif
     endif
    endif
   endif
  endif
 endif
endif


