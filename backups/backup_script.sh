#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

. "$DIR""/../.envs/.production/.django"
. "$DIR""/../.envs/.production/.postgres"

CELERY_BROKER_URL="${REDIS_URL}"
DATABASE_URL="postgres://${POSTGRES_USER}:${POSTGRES_PASSWORD}@${POSTGRES_HOST}:${POSTGRES_PORT}/${POSTGRES_DB}"
LOCAL_BACKUP_DIR="$DIR""/archives/"
REMOTE_BACKUP_DIR="/tmp/"
DATE=`eval date +%Y%m%d`

if docker-compose -f "$DIR""/../production.yml" exec -e DATABASE_URL=$DATABASE_URL -e CELERY_BROKER_URL=$CELERY_BROKER_URL django python manage.py check_studies > /dev/null 2>&1; then
   echo "Writing Postgres dump"
   docker-compose -f "$DIR""/../production.yml" exec postgres sh -c "pg_dump --clean -U $POSTGRES_USER $POSTGRES_DB > /backups/DB_backup"
   mv "$DIR""/DB_backup" "$DIR""/../sdap/media/"
   echo "Building archive"
   tar -pzcf  "$LOCAL_BACKUP_DIR""archive-$DATE"".tar.gz" "$DIR""/../sdap/media/"
   rm "$DIR""/../sdap/media/DB_backup"
fi


cd $LOCAL_BACKUP_DIR
ls -1t  | tail -n +4 | xargs rm -f
echo "Rsync.."
rsync -av --delete --exclude=".*" $LOCAL_BACKUP_DIR $REMOTE_BACKUP_DIR
