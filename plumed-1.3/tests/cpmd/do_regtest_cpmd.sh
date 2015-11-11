#!/bin/bash

###################################################################################
# do_regtest some info [Joost VandeVondele, 2005-02-25 CP2K project]
# [ Davide Branduardi readapted a bit for metadynamics test later on 2008-07-17 ]
#
# Regression testing  ( courtesy of Joost of cp2k team ) 
#    - maintain code quality
#    - helps rapid developemnent and refactoring
#
# What does a regtest do
#    - executes a list of tests
#    - compares the results (outputs) with those of the last known result (reference)
#    - produces a summary
#
# How to set up a regtest
#    - you must be able to build and run on the given machine, the rest should be 'easy'
#    - decide on a directory for doing the regtest, there will be plenty of files in this dir
#      (after a while) so make it something like $HOME/rt
#    - cp $HOME/.....do_regtest $HOME/rt
#    - modify the do_regtest script to match your local environment 
#    - execute './do_regtest' regularly (i.e. after comitting new code)
#
# Interpretation of the results
#  test can be:
#    - 'OK' if the results match those of a previous run precisely. The execution time is also given.
#    - 'NEW' if they have not been executed previously. The reference result is generated 
#      automatically in this run. Tests can also be 'NEW' if they have been reset, i.e. been newly 
#      added to the TEST_FILES_RESET files.
#    - 'RUNTIME FAILURE' if they stopped unexpectedly (e.g. core dump, or stop)
#    - 'WRONG RESULT' if they produce a result that deviates (even a tiny bit) from an old reference
#  the last two options generally mean that a bug has been introduced, which requires investigation.
#  since regtesting only yields information relative to a previously known result, it is most useful
#  to do a regtest before and after you make changes.
#
# Adding/resetting/creating tests to the testsuite
#  these is fully controlled by the following files in the tests directories
#  -TEST_DIRS  : is just a list of directories that contain tests. You can add your directory here( # for comments are allowed ).
#  -TEST_TYPES : this file allows you to create a new test type. I.e. to specify for which words should
#                be grepped and what field should be used in the numerical comparison. 
#                it's format is in the following
#     
#             grep the     del  num      std_out     additional 
#           last line      imi  ber        or        script 
#           containing     ter   of      fileame     to be exec 
#          this keyword        field                 in the dir
#                 !         !    !          !          !
#                 V         V    V          V          V
#
#              ENERGY:      !    6    !  std_out     
#                           !    0    !  HILLS    
#                           !    0    !  HILLS   ! copy_saved_into_restart.sh 
#                
#  -TEST_FILES : the list of input files that need to be executed. You can add your file name here.
#                adding a comment about what it tests might help later debugging problems if a regtest
#                fails
#
#  -TEST_FILES_RESET : you can add files for which the reference output became invalid (e.g. bug fix)
#                      to this list fo files. However be absolutely sure that the change is due to
#                      a bug fix, do not reset these that fail because of unclear reasons. Try to add
#                      a comment to the cvs message and/or the file itself
#
# Command line switches to the do_regtest script (also configurable from within the script)
#  -noreset : do not reset the reference outputs automatically
#  -skipdir string : this switch can repeat, exclude certain dirs from regtesting, useful to 
#                    speed-up regtesting after very localised changes (e.g. -skipdir QS/regtest)
#  -restrictdir string : this switch can repeat, restrict regtesting to certain dirs, useful to 
#                        speed-up regtesting after very localised changes (e.g. -restrictdir QS/regtest)
#
# Script configuration. The value of the follow variables can be redefined, see below
#    awk, datum_full, datum_short,
#    noreset, ndirtoskip, skip_dirs, restrict_dir, ndirtorestrict
#
########################################################################################################
#
# THESE VARIABLES WILL NEED CHANGING FOR YOUR LOCAL CONFIGURATION
#
# - dir_base: the base directory for testing (e.g. /Users/max/build/PLUMED/md_meta/tests/gromacs)
# - maxtasks: how many instances   should run simultaneously (~> #CPUs)
#

dir_base="/home/fabio/data/PLUMED/plumed-cvs-28oct11/md_meta/tests/cpmd"
maxtasks=1
cpmd_prefix="/home/fabio/data/CPMD/CPMD-3.15.1-plumed30sep11/cpmd.x"

#
# The following variables typically need no changes on Linux machine, but might need changes on
# other an OS
#

# *** how to execute an input file 

# *** make and awk
awk=awk
#awk=nawk

# *** a short and long version of the data, in a format that CVS understands
#datum_full=`date --iso-8601="seconds"`
#datum_short=`date --iso-8601="seconds"`
datum_full=`date '+%Y-%m-%dT%H:%M:%S+0100'`
datum_short=`date '+%Y-%m-%d'`

# *** default settings for command line switches
noreset="reset"
ndirstoskip=0
skip_dirs[1]=""
ndirstorestrict=0
restrict_dirs[1]=""

###################################################################################
#
# From here on no changes to the script should be needed
#
###################################################################################
#
# command line argument passing
#
###################################################################################
while [ $# -ge 1 ]; do
case $1 in
  # build a list of directories to skip
  -skipdir) let ndirstoskip=ndirstoskip+1;
        skip_dirs[ndirstoskip]=$2;
        shift;
        ;;
  # build a list of directories to restrict, i.e. only matching dirs will be run 
  -restrictdir) let ndirstorestrict=ndirstorestrict+1;
        restrict_dirs[ndirstorestrict]=$2;
        shift;
        ;;
  # do not reset reference outputs
  -noreset) noreset="noreset";;
  *)  break;;
esac
shift
done

###################################################################################
#
# set up the initial directory structures
#
###################################################################################
test_types_file=${dir_base}/TEST_TYPES
dir_last=${dir_base}/LAST-cpmd
dir_out=${dir_base}/TEST-cpmd
changelog_diff=${dir_out}/ChangeLog.diff
changelog_diff_tests=${dir_out}/ChangeLog-tests.diff
error_description_file=${dir_out}/error_summary
mkdir -p ${dir_out}
mkdir -p ${dir_last}
rm -fR ${error_description_file}
touch  ${error_description_file}
touch  ${dir_last}/ChangeLog-tests

###################################################################################
#
# simple function to end the tests all in the same way
#
###################################################################################
function end_test() {
echo "--------------------------------------------------------------------------"
date
echo "*************************** testing ended ********************************"
exit $1
}

###################################################################################
#
# function to grep for changes in the output. Takes five arguments
#
###################################################################################
function do_test_grep(){
 output_new=$1
 output_old=$2
 error_file=$3
 grep_string=$4
 grep_field=$5 
 grep_file=$5 
 if [ "${grep_field}" == "0" ] ; then  #compare all the line
    grep_string="ALL THE LINE "
    e1=` tail -1 ${output_old}  ` 
    e2=` tail -1 ${output_new}  ` 
    if [ "$e1" == "$e2" ] ; then 
       big=0 
    else
       big="DIFFERENT LINE " 
    fi
 elif [ "${grep_field}" -gt "0"  ] ; then
    e1=`grep -a "${grep_string}" ${output_old} | tail -1 | ${awk} -v f=${grep_field} '{print $f}'`
    e2=`grep -a "${grep_string}" ${output_new} | tail -1 | ${awk} -v f=${grep_field} '{print $f}'`
    big=`echo "${e1} ${e2}" | ${awk} '{v=sqrt((($1-$2)/$2)^2); if (v>1.0E-14) printf("%16.8e",v); else printf("0") ;}'`
 else 
    e1=`grep -A1 -a "${grep_string}" ${output_old} | tail -1 | ${awk} -v f=${grep_field} '{f=-f;print $f}'`
    e2=`grep -A1 -a "${grep_string}" ${output_new} | tail -1 | ${awk} -v f=${grep_field} '{f=-f;print $f}'`
    big=`echo "${e1} ${e2}" | ${awk} '{v=sqrt((($1-$2)/$2)^2); if (v>1.0E-14) printf("%16.8e",v); else printf("0") ;}'`
 fi
 case ${big} in
 0)
  # ok, same energy
  return 0 ;;
 *)
  # nope too large
  echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_file}
  echo "${output_new} : " >> ${error_file}
  echo " ${grep_string} : old = ${e1} new = ${e2} " >> ${error_file}
  echo " relative error : ${big}  " >> ${error_file}
  echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_file}
  return 1 ;;
 esac 
}

###################################################################################
#
# function to select which test to run
#
###################################################################################
function do_test() {
 which_test=$1
 output_new=$2
 output_old=$3
 error_file=$4
 case ${which_test} in
 0) 
   #just be happy you executed
   return 0;;
 *)
   do_test_grep ${output_new} ${output_old} ${error_file} "${test_grep[which_test]}" "${test_col[which_test]}" "${test_file[which_test]}" 
   return $? ;;
 esac
}

# *** start testing
echo "*************************** testing started ******************************"
echo " started on " `date`
echo "-------------------------regtesting cpmd----------------------------------"

###################################################################################
#
# parse the TEST_TYPES file to do different kinds of test 
#
# tests grep for the last line in the file where a string matches (test_grep) 
# and compares a numeric field at a given column (test_col)
#
# the format of the TEST_TYPES file is (notice the '!' as a field separator, to allow
# for spaces in the test_grep)
#
# Ntest_types
# test_grep_1 ! test_col_1
# test_grep_2 ! test_col_2
# ....
# followed by comment lines
#
###################################################################################
Ntest_types=`awk -v l=1 -v c=1 'BEGIN{FS="!"}{lr=lr+1;if (lr==l) print $c}' ${test_types_file}`
test_grep[0]=""
test_col[0]=1
t=1
while [ $t -le ${Ntest_types} ]; do
test_grep[t]=`grep -v "#" ${test_types_file} | ${awk} -v l=$t -v c=1 'BEGIN{FS="!"}{lr=lr+1;if (lr==l+1) printf("%s",$c)}' `
test_col[t]=`grep -v "#"  ${test_types_file} | ${awk} -v l=$t -v c=2 'BEGIN{FS="!"}{lr=lr+1;if (lr==l+1) printf("%d",$c)}' `
test_file[t]=`grep -v "#"  ${test_types_file} | ${awk} -v l=$t -v c=3 'BEGIN{FS="!"}{lr=lr+1;if (lr==l+1) printf("%s",$c)}' | sed s/" "//g `

number_of_fields=`grep -v "#"  ${test_types_file} | ${awk} -v l=$t  'BEGIN{FS="!"}{lr=lr+1;if (lr==l+1) printf("%d",NF)}'  `
# additional script ? 
if [ "${number_of_fields}" -eq "4" ] ; then 
  test_script[t]=`grep -v "#"  ${test_types_file} | ${awk} -v l=$t -v c=4 'BEGIN{FS="!"}{lr=lr+1;if (lr==l+1) printf("%s",$c)}' `
else
  test_script[t]="X";
fi
#echo "TEST TYPES ${test_grep[t]} ${test_col[t]} ${test_file[t]} ${test_script[t]} ${number_of_fields}" 
let t=t+1
done
#exit
###################################################################################
#
# *** now start testing 
# *** for a given directory we do a run on all files in TEST_FILES and
# *** do the test as indicated by the number
# *** files are run in order so that they can e.g. restart
#
###################################################################################
n_runtime_error=0
n_wrong_results=0
n_correct=0
n_tests=0
n_new=0

#
# get a list of directories to be tested, taking care of the exclusions
#
dirs=`cat ${dir_base}/TEST_DIRS | grep -v "#"`
newdirs=""
for dir in ${dirs}
do
  match="no"
  t=1
  # *** match to exclusion list
  while [ $t -le ${ndirstoskip} ]; do
     if [[ "${skip_dirs[t]}" == "${dir}" ]]; then
        match="yes" 
     fi
     let t=t+1
  done
  # *** match to the restrict list, if no restrict list is found, all dirs match
  if [ ${ndirstorestrict} -gt 0 ]; then
     restrictmatch="no"
     t=1
     while [ $t -le ${ndirstorestrict} ]; do
        if [[ "${restrict_dirs[t]}" == "${dir}" ]]; then
           restrictmatch="yes" 
        fi
        let t=t+1
     done
  else
    restrictmatch="yes"
  fi

  # *** if not excluded add to list of dirs
  if [[ "${match}" == "no" && "${restrictmatch}" == "yes" ]]; then
     new_dirs="$new_dirs $dir"
  fi
done
dirs=$new_dirs

#
# execute all regtests
#
cd ${dir_base}
rm -f REGTEST_RUNNING-* REGTEST_TASK_RESULT-* REGTEST_TASK_TESTS-*

for dir in ${dirs};
do
 #
 # tests in different dirs can run in parallel. We spawn processes up to a given maximum
 #
 task=${dir//\//-}
 (
  touch ${dir_base}/REGTEST_RUNNING-$task
  n_runtime_error=0
  n_wrong_results=0
  n_correct=0
  n_tests=0
  n_new=0

  # enter where the tests are located
  cd ${dir_base}/${dir}
  # make the directory in the target location 
  mkdir -p ${dir_out}/${dir}
  # make the directory in the target location for the last reliable result 
  mkdir -p ${dir_last}/${dir}
  touch ${dir_last}/${dir}/TEST_FILES_RESET

  # 
  # first reset reference outputs that have become out-dated since the last run
  # these files must be recalculated brand new
  #
  if [[ ${noreset} != "noreset" ]]; then
     diff TEST_FILES_RESET ${dir_last}/${dir}/TEST_FILES_RESET > ${dir_out}/${dir}/TEST_FILES_RESET.diff
     cp TEST_FILES_RESET ${dir_last}/${dir}/TEST_FILES_RESET
     nreset=`grep '<' ${dir_out}/${dir}/TEST_FILES_RESET.diff | grep -v '#' |  ${awk} '{c=c+1}END{print c}'`
     for ((itest=1;itest<=nreset;itest++));
     do
        reset_file=`grep '<' ${dir_out}/${dir}/TEST_FILES_RESET.diff | grep -v '#' | ${awk} -v itest=$itest '{c=c+1;if (c==itest) print $2}'`
        rm -f ${dir_last}/${dir}/${reset_file}.log
     done
  fi
  #
  # run the tests now
  #
  echo "Starting tests in ${dir_base}/${dir}"
  echo ">>>>>>>>>>>>>>>>> ${dir_base}/${dir}" > ${dir_base}/REGTEST_TASK_TESTS-$task
  if [ ! -e "TEST_FILES" ] ; then
       echo "MISSING TEST_FILES IN ${dir_base}/${dir} : CANNOT PERFORM THE TEST"
       exit
  fi
  ntest=`grep -v "#" TEST_FILES | ${awk} '{c=c+1}END{print c}'`
  for ((itest=1;itest<=ntest;itest++));
  do
     n_tests=$((n_tests+1))
     this_test=""
     input_file=`grep -v "#" TEST_FILES | ${awk} -v itest=$itest '{c=c+1;if (c==itest) print $1}'`
     # just one test right now, but this should generalize
     test_types=`grep -v "#" TEST_FILES | ${awk} -v itest=$itest '{c=c+1;if (c==itest) print $2}'`
     output_file=${dir_out}/${dir}/${input_file}.log
     output_last=${dir_last}/${dir}/${input_file}.log

     # compare file: is the file through which you should compare

     if  [ "${test_file[test_types]}" == "std_out"  ] ; then 
         compare_file=${dir_out}/${dir}/${input_file}.log
         compare_last=${dir_last}/${dir}/${input_file}.log
     else
         compare_file=${dir_out}/${dir}/${test_file[test_types]}.log
         compare_last=${dir_last}/${dir}/${test_file[test_types]}.log #to conformate with reset syntax 
     fi
     #echo "MY TESTFILE IS  ${test_file[test_types]} FILE ${compare_file} LAST ${compare_last} " 

     #  
     # put everything in the output file
     #  
     if [ "${test_script[test_types]}" != "X" ] ; then
           echo " executing ${test_script[test_types]} first "
           ${test_script[test_types]}
     fi 

     ($cpmd_prefix ${input_file} . -plumed > ${output_file} 2>/dev/null )

     #  
     # in case you don't have to compare the output
     #  
     if  [ "${test_file[test_types]}" != "std_out"  ] ; then 
        cp ${test_file[test_types]}  ${compare_file} 
     fi
     # *** cpmd failed obviously
     if (( $? )); then
        echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
        echo ${output_file} >> ${error_description_file}
        tail -40 ${output_file} >> ${error_description_file}
        echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
        this_test="RUNTIME FAIL"
        n_runtime_error=$((n_runtime_error+1))
        failed_tests="${failed_tests} ${output_file}"
     else 
        # *** but didn't end !? (namd2)
        grep -a " " ${output_file} &> /dev/null
        if (( $? )); then
           echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
           echo ${output_file} >> ${error_description_file}
           tail -40 ${output_file} >> ${error_description_file}
           echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
           this_test="RUNTIME FAIL"
           n_runtime_error=$((n_runtime_error+1))
           failed_tests="${failed_tests} ${output_file}"
        else
           # *** still running, you must be joking...
           # *** see if we manage to pass the testing
           # *** but only if we can compare
           if [ -f ${compare_last} ]; then
              for test_type in ${test_types};
              do
                 do_test ${test_type} ${compare_file} ${compare_last} ${error_description_file}
                 if (( $? )); then
                    this_test="WRONG RESULT TEST ${test_type}"
                    n_wrong_results=$((n_wrong_results+1))
                    # *** no further testing
                    break;
                 else
                    n_correct=$((n_correct+1))
                    this_test="OK"
                 fi
              done
           else
              this_test="NEW"
              n_new=$((n_new+1))
           fi
        fi
     fi
     # Keep the output up-to-date
 
     case ${this_test} in
     "NEW" )
        #echo "NOW COPY ${compare_file} INTO  ${compare_last} "
        cp ${compare_file} ${compare_last} 
        # cpmd specific timing 
        timing=`grep "TOTAL TIME" ${output_file} | ${awk} '{printf("%6.2f",$3)}'`
        this_test="${this_test}  (${timing} sec)" ;;
     "OK" )
        # cpmd specific timing 
        timing=`grep "TOTAL TIME" ${output_file} | ${awk} '{printf("%6.2f",$3)}'`
        this_test="${this_test}  (${timing} sec)" ;;
         
     esac
     type=" TYPE ${test_file[${test_type}]} "
 
     printf "%50s  %20s %20s\n" "${dir}/${input_file}"  "${this_test}" "${type}" >> ${dir_base}/REGTEST_TASK_TESTS-$task
  done
  echo "<<<<<<<<<<<<<<<<< ${dir_base}/${dir}" >> ${dir_base}/REGTEST_TASK_TESTS-$task
  echo "${n_runtime_error} ${n_wrong_results} ${n_correct} ${n_new} ${n_tests}" > ${dir_base}/REGTEST_TASK_RESULT-$task
  cat ${dir_base}/REGTEST_TASK_TESTS-$task
  rm -f ${dir_base}/REGTEST_TASK_TESTS-$task ${dir_base}/REGTEST_RUNNING-$task
 )&

 #
 # here we allow only a given maximum of tasks
 #
 runningtasks=10000
 while (( runningtasks >= maxtasks ))
 do
   sleep 1
   runningtasks=`ls -1 ${dir_base}/REGTEST_RUNNING-* 2> /dev/null | awk 'BEGIN{c=0}{c=c+1}END{print c}'`
 done

done

#
# wait for all tasks to finish
#
wait
#
# generate results
#
for dir in ${dirs};
do
  task=${dir//\//-}
  file=${dir_base}/REGTEST_TASK_RESULT-$task
  tmp=`awk '{print $1}' $file`
  n_runtime_error=$((n_runtime_error+tmp))
  tmp=`awk '{print $2}' $file`
  n_wrong_results=$((n_wrong_results+tmp))
  tmp=`awk '{print $3}' $file`
  n_correct=$((n_correct+tmp))
  tmp=`awk '{print $4}' $file`
  n_new=$((n_new+tmp))
  tmp=`awk '{print $5}' $file`
  n_tests=$((n_tests+tmp))
  rm -f $file
done

echo "--------------------------------------------------------------------------"
cat "${error_description_file}"
echo "--------------------------------- summary --------------------------------"
printf "number of FAILED  tests %d\n" ${n_runtime_error}
printf "number of WRONG   tests %d\n" ${n_wrong_results}
printf "number of CORRECT tests %d\n" ${n_correct}
printf "number of NEW     tests %d\n" ${n_new}
printf "number of         tests %d\n" ${n_tests}


end_test 0
