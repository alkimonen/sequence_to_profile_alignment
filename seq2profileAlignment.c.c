#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#define MAX_ALINGMENT_LENGTH 512

static struct option long_options[] = {
    { "aln", required_argument, 0, 'a'},
    { "fasta", required_argument, 0, 's'},
    { "gap", required_argument, 0, 'g'},
    { "match", required_argument, 0, 'm'},
    { "mismatch", required_argument, 0, 'p'},
    { "out", required_argument, 0, 'o'}
};

enum directions { M, SG, AG};

char keys[5] = { 'A', 'C', 'G', 'T', '-'};

int matchScore = 1;
int gapPen = -1;
int mmPen = -1;

int noOfAln = 0;
char **alnNames;
char **aln;
int alnLen = 0;

char *seqName;
char *seq;
int seqLen;

double **profile;

double **scores;
enum directions **moves;

char **resultNames;
char **results;
int resLen;
int noOfSeqs;

double calcMatchScore( char ch, int index) {
    double score = 0;
    for ( int k = 0; k < 5; k++) {
        if ( ch == keys[k])
            score = score + matchScore * profile[index][k];
        else
            score = score + mmPen * profile[index][k];
    }
    return score;
}

void globalAlignment() {
    scores[0][0] = 0;
    moves[0][0] = M;

    for ( int i = 1; i < alnLen+1; i++) {
        scores[i][0] = scores[i-1][0] + gapPen;
        moves[i][0] = SG;
    }

    for ( int j = 1; j < seqLen+1; j++) {
        scores[0][j] = scores[0][j-1] + gapPen;
        moves[0][j] = AG;
    }

    for ( int i = 1; i < alnLen+1; i++) {
        for ( int j = 1; j < seqLen+1; j++) {
            double maxScore, currentScore;
            enum directions move;

            maxScore = scores[i-1][j-1] + calcMatchScore( seq[j-1], i-1);
            move = M;

            currentScore = scores[i-1][j] + gapPen;
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = SG;
            }

            currentScore = scores[i][j-1] + gapPen;
            if ( currentScore > maxScore) {
                maxScore = currentScore;
                move = AG;
            }

            scores[i][j] = maxScore;
            moves[i][j] = move;
        }
    }
}

void backtraceAlignment() {
    int i, j;
    noOfSeqs = noOfAln + 1;

    resLen = 0;
    i = alnLen;
    j = seqLen;
    while ( i >= 0 && j >= 0) {
        if ( moves[i][j] == M) {
            i--;
            j--;
        }
        else if ( moves[i][j] == SG) {
            i--;
        }
        else if ( moves[i][j] == AG) {
            j--;
        }
        resLen++;
    }
    resLen--;

    results = (char **) malloc( noOfSeqs * sizeof(char *));
    for ( int i = 0; i < noOfSeqs; i++) {
        results[i] = (char *) malloc( resLen);
    }

    i = alnLen;
    j = seqLen;
    for ( int index = resLen-1; index >= 0; index--) {
        if ( moves[i][j] == M) {
            i--;
            j--;
            for ( int n = 0; n < noOfAln; n++) {
                results[n][index] = aln[n][i];
            }
            results[noOfSeqs-1][index] = seq[j];
        }
        else if ( moves[i][j] == SG) {
            i--;
            for ( int n = 0; n < noOfAln; n++) {
                results[n][index] = aln[n][i];
            }
            results[noOfSeqs-1][index] = '-';
        }
        else if ( moves[i][j] == AG) {
            j--;
            for ( int n = 0; n < noOfAln; n++) {
                results[n][index] = '-';
            }
            results[noOfSeqs-1][index] = seq[j];
        }
    }
}

long ms( struct timeval t) {
    return t.tv_sec * 1000000 + t.tv_usec;
}

int main( int argc, char** argv)
{
	char *alnFile, *seqFile, *outFile;
    int index;
    int opt;
    int aFlag = 0, sFlag = 0, oFlag = 0;


	while ((opt = getopt_long( argc, argv, "a:s:g:m:p:o", long_options, &index)) != -1) {
        switch(opt) {
            case 'a':
            	alnFile = (char *) malloc( 25);
            	strcpy( alnFile, optarg);
            	aFlag = 1;
                break;
            case 's':
            	seqFile = (char *) malloc( 25);
            	strcpy( seqFile, optarg);
            	sFlag = 1;
                break;
            case 'g':
                gapPen = atoi( optarg);
                break;
            case 'm':
                matchScore = atoi( optarg);
                break;
            case 'p':
                mmPen = atoi( optarg);
                break;
            case 'o':
            	outFile = (char *) malloc( 25);
            	strcpy( outFile, optarg);
            	oFlag = 1;
                break;
        }
    }

    if ( aFlag == 0 || sFlag == 0 || oFlag == 0) {
        printf( "Arguments missing.\n");
        return 1;
    }

    // allocate memory
    alnNames = (char **) malloc( 10 * sizeof(char *));
    aln = (char **) malloc( 10 * sizeof(char *));

    // read alignment file
    char *word = (char *) malloc( MAX_ALINGMENT_LENGTH);

    FILE *fp = fopen( alnFile, "r");
    if ( fp != NULL) {
        while ( fscanf( fp, "%s", word) != EOF) {
            int len = 0;
            for ( int j = 0; j < MAX_ALINGMENT_LENGTH; j++) {
                if ( !word[j])
                    break;
                len++;
            }
            alnNames[noOfAln] = (char *) malloc( len+1);
            strcpy( alnNames[noOfAln], word);

            if ( fscanf( fp, "%s", word) != EOF) {
                if ( alnLen == 0) {
                    alnLen = 0;
                    for ( int j = 0; j < MAX_ALINGMENT_LENGTH; j++) {
                        if ( !word[j])
                            break;
                        alnLen++;
                    }
                }
                aln[noOfAln] = (char *) malloc( alnLen+1);
                strcpy( aln[noOfAln], word);
            }
            noOfAln++;
        }
    }
    else {
        printf( "Can't open alignment file\n");
        return 1;
    }
    fclose( fp);
    free( alnFile);

    // read sequence file
    char *line;
    size_t len = 0;
    fp = fopen( seqFile, "r");
    if ( fp == NULL) {
        printf( "Can't open sequence file.\n");
        return 1;
    }
    while ( fscanf( fp, "%s", word) != EOF) {
        if ( word[0] == 62) {
            int nameLen = 0;
            for ( int j = 0; j < MAX_ALINGMENT_LENGTH; j++) {
                if ( !word[j])
                    break;
                nameLen++;
            }
            seqName = (char *) malloc( nameLen+1);
            for ( int i = 0; i < nameLen-1; i++)
            {
                seqName[i] = word[i+1];
            }
        }
        else {
            seqLen = 0;
            for ( int j = 0; j < MAX_ALINGMENT_LENGTH; j++) {
                if ( !word[j])
                    break;
                seqLen++;
            }
            seq = (char *) malloc( seqLen+1);
            strcpy( seq, word);
        }
    }
    fclose( fp);
    free( word);
    free( seqFile);

    // create profile
    profile = (double **) malloc( alnLen * sizeof( double *));
    for ( int i = 0; i < alnLen; i++) {
        profile[i] = (double *) malloc( 5 * sizeof( double));

        for ( int k = 0; k < 5; k++) {
            double chCount = 0;

            for ( int j = 0; j < noOfAln; j++) {
                if ( aln[j][i] == keys[k])
                    chCount++;
            }
            profile[i][k] = chCount / noOfAln;
        }
    }

    // calculate scores
    scores = (double **) malloc( (alnLen+1) * sizeof( double *));
    moves = (enum directions **) malloc( (alnLen+1) * sizeof(enum directions *));
    for ( int i = 0; i < alnLen+1; i++) {
        scores[i] = (double *) malloc( (seqLen+1) * sizeof( double));
        moves[i] = (enum directions *) malloc( (seqLen+1) * sizeof( enum directions));
    }

    globalAlignment();

    // backstrace alignemnt
    backtraceAlignment();

    // write results into file
    fp = fopen( outFile, "w");
    if ( fp == NULL) {
        printf( "Unable to open output file.\n");
        return 1;
    }

    resultNames = (char **) malloc( noOfSeqs * sizeof(char *));

    for ( int n = 0; n < noOfAln; n++) {
        resultNames[n] = alnNames[n];
    }
    resultNames[noOfSeqs-1] = seqName;

    for ( int n = 0; n < noOfSeqs; n++) {
        fprintf( fp, "%-10s %s\n", resultNames[n], results[n]);
    }

    fclose( fp);
    free( outFile);

    // deallocate memory
    for ( int i = 0; i < noOfAln; i++) {
        if ( alnNames[i])
            free( alnNames[i]);
        if ( aln[i])
            free( aln[i]);
    }
    free( alnNames);
    free( aln);

    free( seqName);
    free( seq);

    for ( int i = 0; i < alnLen; i++) {
        if ( profile[i])
            free( profile[i]);
    }
    free( profile);

    for ( int i = 0; i < alnLen+1; i++) {
        if ( scores[i])
            free( scores[i]);
        if ( moves[i])
            free( moves[i]);
    }
    free( scores);
    free( moves);

    for ( int i = 0; i < noOfSeqs; i++) {
        if ( results[i])
            free( results[i]);
    }
    free( results);
    free( resultNames);

	printf( "Alignment completed.\n");
    return 0;
}
