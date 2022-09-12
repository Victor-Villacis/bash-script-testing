#!/bin/bash

echo "Hello"

# Run date and hostname command and store output to shell variables
now="$(date)"
computer_name="$(hostname)"

#
# print it or use the variable
# Variable names are case sensitive $now and $NOW are different names
#
echo "$var"
echo "Current date and time : $now"
echo "Computer name : $computer_name"
echo ""

# Read input using read command
read -p "May I know your name please? " name
echo "Hello $name, let us be friends"
echo ""
