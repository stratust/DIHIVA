#!/usr/bin/env perl

use strict;
use warnings;
use MyApp;
use Log::Any::App '$log', -screen => 1;    # turn off screen logging explicitly


MyApp->new_with_command->run();
