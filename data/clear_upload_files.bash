file_path="/data/groups/ghub/tools/gisplot2"
if [ "$(ls -A $file_path)" ]; then
    for session in $file_path/*
    do
        #echo $session
        mod=$(date -r $session +%s)
        now=$(date +%s)
        days=$(expr \( $now - $mod \) / 86400)
        #echo $days
        max_exist=3
        if [ $days -gt $max_exist ]; then
            rm -rf $session
        fi
    done
fi