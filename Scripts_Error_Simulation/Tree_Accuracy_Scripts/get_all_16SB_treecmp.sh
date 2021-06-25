#!/bin/bash

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0-0.1/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0-0.1/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0-0.1

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.1-0.2/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.1-0.2/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.1-0.2

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.2-0.3/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.2-0.3/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.2-0.3

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.3-0.4/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.3-0.4/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.3-0.4

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.4-0.5/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.4-0.5/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.4-0.5

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.5-0.6/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.5-0.6/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.5-0.6

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.6-0.7/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.6-0.7/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.6-0.7

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.7-0.8/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.7-0.8/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.7-0.8


rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.8-0.9/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.8-0.9/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.8-0.9

rm Original_Trees/*
rm Error_Trees/*
cp ../../16S.B_Original/TREES_0.9-1/* Original_Trees
cp ../../16S.B_Errors/ERROR_TREES_0.9-1/* Error_Trees
./getTreeCmpData.sh Error_Trees
mv FORMAT_TREECMP_OUTPUT TREECMP_OUTPUT_0.9-1

