#!/bin/bash

# List of observatories
observatories=("loao" "dnsm" "cbnuo" "khao" "maao" "doao" "SAO_C5A150M" "kct_stx16803" "rasa36" "lsgt_asi1600mm") # Add more observatory names as needed

# Start tmux session, reuse if it already exists
tmux has-session -t gppy 2>/dev/null

if [ $? != 0 ]; then
  echo "Creating a new tmux session 'astronomy'..."
  tmux new-session -d -s gppy -n setup # Create a first window (temporary) for setup
else
  echo "tmux session 'gppy' already exists. Attaching..."
  tmux attach-session -t gppy
  exit 0
fi

# Create tmux windows for each observatory and execute commands
for obs in "${observatories[@]}"; do
  echo "Setting up window for $obs..."
  tmux new-window -t gppy -n "$obs"
  tmux send-keys -t "$obs" "ulimit -u 65535" C-m
  tmux send-keys -t "$obs" "conda activate gppy" C-m
  tmux send-keys -t "$obs" "cd ~/gppy" C-m
  tmux send-keys -t "$obs" "python gpwatch.py $obs 2 &" C-m # Run the Python script in the background
done

# Remove the first temporary window
tmux select-window -t gppy:setup
tmux kill-window -t setup

# Attach to the tmux session after all windows have been set up
tmux attach-session -t gppy

