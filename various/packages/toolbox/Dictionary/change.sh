#!/bin/bash
for file in ./*
do
    if test -f $file
    then
        echo $file is document && chmod 664 $file
    fi
    if test -d $file
    then
        echo $file is dirctory && chmod 775 $file && cd $file
        for file in ./*
        do
            if test -f $file
            then 
                echo $file is document && chmod 664 $file
            fi
            if test -d $file
            then 
                chmod 775 $file && cd $file
                for file in ./*
                do
                    if test -f $file
                    then 
                        chmod 664 $file
                    fi
                    if test -d $file
                    then
                        chmod 775 $file && cd $file
                        for file in ./*
                        do 
                            if test -f $file
                            then
                                chmod 664 $file
                            fi
                            if test -d $file
                            then
                                chmod 775 $file
                            fi
                        done 
                        cd ..
                    fi
                done
                cd ..
            fi
        done
        cd ..
    fi
done

chmod 775 change.sh
