#!/bin/sh
D1CLASSES="@cfg_branch @cfg_choice @cfg_const @cfg_entry @cfg_files @cfg_menu @cfg_repeat"
D2CLASSES="@cfg_exbranch"
# @cfg_branch is the master in-tree class
INCLASSES="@cfg_choice @cfg_repeat"

echo "Copying @cfg_item/subsasgn.m to all derived classes"
for DC in $D1CLASSES $D2CLASSES; do
    cp @cfg_item/subsasgn.m $DC;
done

echo "Copying @cfg_item/subs_fields.m to all derived classes"
for DC in $D1CLASSES $D2CLASSES; do
    cp @cfg_item/subs_fields.m $DC;
done

echo "Copying @cfg_item/subsref.m to all derived classes"
for DC in $D1CLASSES $D2CLASSES; do
    cp @cfg_item/subsref.m $DC;
done

echo "Copying @cfg_branch/all_leafs.m to all in-tree classes"
for IC in $INCLASSES; do
    cp @cfg_branch/all_leafs.m $IC;
done

echo "Copying @cfg_branch/all_set.m to all in-tree classes"
for IC in $INCLASSES; do
    cp @cfg_branch/all_set.m $IC;
done

echo "Copying @cfg_branch/list.m to all in-tree classes"
for IC in $INCLASSES; do
    cp @cfg_branch/list.m $IC;
done

echo "Copying @cfg_branch/tagnames.m to all in-tree classes"
for IC in $INCLASSES; do
    cp @cfg_branch/tagnames.m $IC;
done

echo "Copying @cfg_choice/clearval.m to @cfg_repeat"
cp @cfg_choice/clearval.m @cfg_repeat/;