#!/bin/bash
# 重写历史提交的 author/committer 邮箱
export FILTER_BRANCH_SQUELCH_WARNING=1
OLD_EMAIL="300889223+Csu-MIB@users.noreply.github.com"
CORRECT_NAME="zzc152"
CORRECT_EMAIL="126441996+zzc152@users.noreply.github.com"

git filter-branch -f --env-filter '
if [ "$GIT_COMMITTER_EMAIL" = "'"$OLD_EMAIL"'" ]; then
    export GIT_COMMITTER_NAME="'"$CORRECT_NAME"'"
    export GIT_COMMITTER_EMAIL="'"$CORRECT_EMAIL"'"
fi
if [ "$GIT_AUTHOR_EMAIL" = "'"$OLD_EMAIL"'" ]; then
    export GIT_AUTHOR_NAME="'"$CORRECT_NAME"'"
    export GIT_AUTHOR_EMAIL="'"$CORRECT_EMAIL"'"
fi
' --tag-name-filter cat -- --branches --tags
