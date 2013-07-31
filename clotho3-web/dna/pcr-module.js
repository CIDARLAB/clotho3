'use strict';

Application.Dna.service('PCR', ['Clotho', 'DNA', 'restrictionSite', function(Clotho, DNA, restrictionSite) {

    /**************
     Primer Design
     **************/

    /**
     * @description Get clamp primer using beginning of sequence
     * @param sequence
     * @param minTemp Temperature (°C) primer should be. Default 50°C
     * @returns {string} primer reverse complimentary to sequence, maximum length 50nt
     */
    var generatePrimer = function (sequence, minTemp) {
        minTemp = !!minTemp ? minTemp : 50;

        var index = minTemp / 3.5, //assumes GC_content < 75%
            primer = sequence.substring(0, index);

        while( DNA.melting_temp_basic(primer) < minTemp || index > 50){
            index = index +1;
            primer = sequence.substring(0, index);
        }

        return DNA.revcomp(primer);
    };

    /**************
     Annealing
     **************/



    var primerBinding = function (primer, sequence) {
        //todo - move to fuzzy search


    };



    /**************
     Alignment
     **************/




    /**************
     Construction
     **************/


    var ligate = function (fragments) {

    };


    return {

    }
}]);