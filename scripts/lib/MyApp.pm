#!/usr/bin/env perl

package MyApp;
use feature qw(say);
use MooseX::App qw(Color);
use Log::Any '$log';
use MooseX::FileAttribute;

app_exclude 'MyApp::Role::JSON';

has 'log' => (
    is            => 'ro',
    isa           => 'Log::Any::Proxy',
    required      => 1,
    default       => sub { Log::Any->get_logger },
    documentation => 'Keep Log::Any::App object',
);

__PACKAGE__->meta->make_immutable;

1;
