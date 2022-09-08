#
# The SubspeciesClassification application
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run qw(run);
use Cwd;
use Clone;

my $script = Bio::KBase::AppService::AppScript->new(\&process_subspeciesclass, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    print STDERR "preflight subspeciesclass ", Dumper($params, $app);

    my $token = $app->token();
    my $ws = $app->workspace();

    # TODO (ask bob): estimate cpu, memory, and runtime values
    # TODO create group of genomes for testing
}

sub process_subspeciesclass
{
    my($app, $app_def, $raw_params, $params) = @_;

    print 'Proc subspecies classification ', Dumper($app_def, $raw_params, $params);

    my $token = $app->token();
    my $ws = $app->workspace();

    #
    # Create an output directory under the current dir. App service is meant to invoke
    # the app script in a working directory; we create a folder here to encapsulate
    # the job output.
    #
    # We also create a staging directory for the input files from the workspace.
    #

    # TODO: may not need a staging directory

    # my $cwd = getcwd();
    my $cwd = File::Temp->newdir( CLEANUP => 1 );
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";

    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => $data_api };
    my $sstring = encode_json($dat);

    # TODO: is this needed
    my $params_to_app = Clone::clone($params);

    #
    # Write job description.
    #
    my $jdesc = "$cwd/jobdesc.json";
    open(JDESC, ">", $jdesc) or die "Cannot write $jdesc: $!";
    print JDESC JSON::XS->new->pretty(1)->encode($params_to_app);
    close(JDESC);

    my $parallel = $ENV{P3_ALLOCATED_CPU};

    my @cmd = ("run_subspecies_classification","-o",$work_dir,"--jfile", $jdesc);

    warn Dumper (\@cmd, $params_to_app);

    my $ok = run(\@cmd);

    if (!$ok)
    {
        die "Command failed: @cmd\n";
    }

    my @output_suffixes = ([qr/\.tsv$/, 'tsv'],[qr/\.json$/, 'json'],[qr/\.tre$/, 'nwk']);

    my $outfile;
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    # TODO: not sure what this does?
    my @files = sort {$a cmp $b } grep { -f "$work_dir/$_" } readdir(D);

    my $output = 1;
    my $output_dir = "$params->{output_path}/.$params->{output_file}";
    for my $file (@files)
    {
        for my $suf (@output_suffixes)
        {
            if ($file =~ $suf->[0])
            {
                $output = 0;
                my $type = $suf->[1];

                $app->workspace->save_file_to_file("$work_dir/$file", {}, "$output_dir/$file", $type, 1,
                                                    (-s "$work_dir/$file" > 10_000 ? 1 : 0), #use shock for larger files
                                                    $token);
            }
        }
    }
}
