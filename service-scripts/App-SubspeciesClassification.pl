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

    return {
	cpu => 2,
	memory => "64G",
	runtime => 10800,
	storage => 0,
    };
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

    my %suffix_map = (tsv => 'tsv',
                      result => 'tsv',
                      err => 'tsv',
                      tre => 'nwk',
                      html => 'html',
                      txt => 'txt',
                      json => 'json');

    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

    if (opendir(my $dh, $work_dir))
    {
        while (my $p = readdir($dh))
        {
            next if $p =~ /^\./;

            my @cmd = ("p3-cp", "-r", "-f", @suffix_map, "$work_dir/$p", "ws:" . $app->result_folder);
            print "@cmd\n";
            my $ok = IPC::Run::run(\@cmd);
            if (!$ok)
            {
                warn "Error $? copying output with @cmd\n";
            }
        } 
        closedir($dh);
    }
    else
    {
        warn "Output directory $work_dir does not exist\n";
    }
}
