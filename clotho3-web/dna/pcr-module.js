"use strict";

Application.Dna.service('PCR', ['DNA', 'Digest', function(DNA, Digest) {

    /*************
     Classes
     ************/

    //note - assumes only terminal ends -- run through parse first
    function Fragment (markedSeq) {
        /* checks */
        //already a fragment - just pull the sequence and reprocess
        if (markedSeq instanceof Fragment) {
            markedSeq = markedSeq.sequence;
        }
        //nonterminal marks --- shouldn't be passed in...
        if ((Digest.findOverhangs(markedSeq, true)).length) {
            //future - how to handle?
        }

        this.sequence = markedSeq;
        this.process();
    }

    Fragment.prototype = {

        /**** prototype ****/

        process : function(newSequence) {
            if (!!newSequence)
                this.sequence = newSequence;
            this.ends = Digest.findOverhangs(this.sequence);
            this.endMatches = _.pluck(this.ends, 'match');
        },
        reverse : function() {
            this.sequence = DNA.revcomp(this.sequence);
            this.process()
        },
        circularize : function() {
            if (this.endMatches.length == 2 && this.endsMatch(this.endMatches[0], this.endMatches[1])) {
                //remove end marks -- remove last overhang to circularize, and remove cut marks in first terminal cut

                //remove last overhang
                this.sequence = this.sequence.substring(0, this.sequence.length - this.endMatches[1].length);
                //remove first cut marks, don't remove in whole sequence
                this.sequence = this.ends[0].overhang + this.sequence.substring(this.endMatches[0].length);
                this.process();
            }
        },

        /**** basic matching ****/

        //check that ends are complementary
        endsMatch : function(match1, match2) {
            return (match1 == match2 || match1 == DNA.revcomp(match2))
        },
        //only checks revcomp match
        endsMatchRevcomp : function (match1, match2) {
            return match1 == DNA.revcomp(match2)
        },

        /************
        matching
        ***********/

        /***** match string ****/

        //given match, returns array of ends in this that match
        findMatchingEnds : function(match) {
            return _.filter(this.ends, function(end, index) {
                return this.endsMatch(end.match, match);
            }, this);
        },
        //given match, returns array of indices of ends in this that match
        findMatchingEndIndices : function(match) {
            var endsMatching = this.findMatchingEnds(match);

            switch (endsMatching.length) {
                case 0 : {
                    return false;
                }
                case 1 : {
                    return [_.indexOf(this.ends, endsMatching[0])]
                }
                default : {
                    //because assume <= 2 ends per fragment
                    return [0,1];
                }
            }
        },
        //given match, whether fragment has a matching end in this fragment
        //can also pass in fragment, defers to canMatchFrag
        canMatch : function (match) {
            if (match instanceof Fragment) {
                return this.canMatchFrag(match)
            } else {
                return !!(this.findMatchingEnds(match)).length;
            }
        },

        /**** match fragment ****/

        canMatchFrag : function (otherFrag) {
            return !!(_.find(otherFrag.endMatches, function (otherMatch) {
                return this.canMatch(otherMatch);
            }, this) || []).length;
        },

        //return object of ends which match: keys = this.ends index, values = otherFrag
        matchMap : function (otherFrag) {
            var connections = {};
            _.each(this.endMatches, function (thisMatch, thisInd) {
                _.each(otherFrag.endMatches, function (otherMatch, otherInd) {
                    if (this.endsMatch(thisMatch, otherMatch)) {
                        if (!connections[thisInd])
                            connections[thisInd] = otherInd;
                        else {
                            //points to two fragments
                            //future - handle this
                            console.log('pointing to two fragments');
                        }
                    }
                }, this)
            }, this);

            if (!_.keys(connections).length) return undefined;
            return connections
        },

        /*** match array of fragments***/

        matchMapArray : function(otherFrags) {

            //don't introduce return inconsistency
            /*
            if (otherFrags.length == 1) {
                return this.matchMap(otherFrags[0])
            } else if (otherFrags instanceof Fragment) {
                return this.matchMap(otherFrags)
            }
            */

            var matchMap = {};
            _.each(otherFrags, function(frag, index) {
                if (this !== frag)
                    matchMap[index] = this.matchMap(frag);
            }, this);
            return matchMap;
        },
        findFirstMatch : function (otherFrags) {
            return _.find(otherFrags, function(otherFrag, index) {
                //this !== otherFrag && console.log(this, otherFrag); //testing
                return ((this !== otherFrag) && (this.canMatch(otherFrag)))
            }, this)
        },
        findFirstMatchIndex : function (otherFrags) {
            return _.indexOf(otherFrags, this.findFirstMatch(otherFrags));
        },
        canMatchArray : function (otherFrags) {
            return !!this.findFirstMatch(otherFrags)
        },

        /*** joining ***/

        joinFragment : function(otherFrag) {
            if (_.isString(otherFrag))
                otherFrag = new Fragment(otherFrag);

            var connections = this.matchMap(otherFrag);

            if (!_.keys(connections).length)
                return false;

            //find connection
            //for now just pull the first connection
            //future - check for unique matches
            var firstConnection = _.pairs(connections)[0],
                thisEndInd = firstConnection[0],
                otherEndInd = firstConnection[1],
                thisEnd = this.ends[thisEndInd],
                otherEnd = otherFrag.ends[otherEndInd];
            
            //console.log('joinFragment - this end and other end:', thisEnd, otherEnd);

            //console.log('before orienting', thisEndInd, otherEndInd);


            //todo - own function, standard return to avoid all these variables
            //orient fragments
            //can assume terminal, so just need to make sure one is zero
            if (thisEnd.index == 0) {
                if (otherEnd.index == 0) {
                    otherEndInd = (!!otherEndInd || otherFrag.ends.length == 1) ? 0 : 1;
                    otherFrag.reverse();
                    otherEnd = otherFrag.ends[otherEndInd];
                }
            } else {
                if (otherEnd.index != 0) {
                    otherEndInd = (!!otherEndInd || otherFrag.ends.length == 1) ? 0 : 1;
                    otherFrag.reverse();
                    otherEnd = otherFrag.ends[otherEndInd];
                }
            }

            //console.log('after orienting', thisEndInd, otherEndInd);

            var firstFrag = (thisEnd.index == 0) ? otherFrag : this,
                firstEnd = (thisEnd.index == 0) ? otherEnd : thisEnd,
                secondFrag = (firstFrag === this) ? otherFrag : this,
                secondEnd = (firstEnd === thisEnd) ? otherEnd : thisEnd;

            //console.log('first frag', firstFrag, 'first end', firstEnd, 'second frag', secondFrag, 'second end', secondEnd);

            //join
            //todo - options for html and cut marks etc.

            var joined = firstFrag.sequence.substring(0, firstEnd.index) +
                Digest.removeMarks(firstEnd.match) +
                secondFrag.sequence.substring(secondEnd.index + secondEnd.match.length);

            //console.log(joined);


            //note - can't leave in cut marks and process -- just use directive tags
            this.process(joined);
            return true
        },
        alignFragment : function (otherFrag) {

        }
    };


    /**************
     Primer Design
     **************/

    /**
     * @description Get clamp primer using beginning of sequence
     * @param sequence
     * @param minTemp Temperature (°C) primer should be. Default 50°C
     * @returns {string} primer reverse complimentary to sequence, minimum 10nt, maximum 50nt
     */
    var generatePrimer = function generatePrimer(sequence, minTemp) {
        minTemp = !!minTemp ? minTemp : 50;

        var minLength = 10,
            maxLength = 50,
            startTest = Math.floor(minTemp / 3.5), //assumes GC_content < 75%
            index = (startTest > minLength ? startTest : minLength),
            primer = sequence.substring(0, index);

        while( DNA.melting_temp_basic(primer) < minTemp || index > maxLength){
            index = index +1;
            primer = sequence.substring(0, index);
        }

        return DNA.revcomp(primer);
    };

    /**************
     Annealing
     **************/

    /**
     * @description Finds indices where primer anneals to sequence, as given by Digest.findIndices() -- i.e. exact match without tail
     * @param sequence
     * @param primer
     * @returns {Object} In form {forward: Array, reverse: Array}} where each array is indices the primer anneals, empty if no matches.
     */
    var findAnnealAllExact = function findAnnealAllExact(sequence, primer) {

        //future - move to fuzzy search, account for tail on primer

        var forward = Digest.findIndices(sequence, primer, false);
        var reverse = Digest.findIndices(sequence, DNA.revcomp(primer), false);

        return {forward: forward, reverse: reverse}
    };


    var findAnnealFuzzy = function(sequence, primer) {
        //todo
    };

    //todo - return array of objects, each with index and if forward and has overhang
    var findAnnealSimple = function(template, primer) {
        //start from 3', go back until have unique match
        var searchFrag, searchReg, result = null;

        for(var initialBack = 8,
                start = primer.length - initialBack,
                matches = {};
            start > 0,
                searchFrag = primer.substring(start),
                searchReg = Digest.createRegex(searchFrag);
            --start)
        {
            console.log(start);

            //check forward
            matches.forward = template.match(searchReg) || [];
            //check revcomp
            matches.reverse = DNA.revcomp(template).match(searchReg) || [];

            console.log(matches, matches.forward.length, matches.reverse.length, matches.forward.length + matches.reverse.length);

            if (!matches.forward.length && !matches.reverse.length) {
                console.log('no *exact* matches found for length ' + start + ' from 3 prime end');
                //todo - step back one, check for single match
                result = false;
                break;
            } else if ((matches.forward.length + matches.reverse.length) == 1) {
                console.log('one match' );
                //todo - go back as far as possible
                //return index of single match
                result = matches.forward.length ? template.search(searchReg) : DNA.revcomp(template).search(searchReg);
                break;
            } else {
                console.log('else');
            }
        }

        return result;
    };

    //wrapper function...
    var anneal = function (template, primer, fuzzy) {
        fuzzy = !!fuzzy || true;

        //future - once write fuzzy, implement here
        return (fuzzy) ? findAnnealSimple(template, primer) : findAnnealSimple(template, primer);

    };

    //e.g. PCA, given array of oligos, find which match which
    var annealArray = function (oligos) {

        var track = {};
        for (var i = 0; i < oligos.length; i++) {
            track[i] = {oligo : oligos[i]};

            //todo - find matches in oligos
        }

    };


    //TODO - write algorithm that runs the PCR, not as predictive

    /**************
     PCR algorithms
     ****************/


    /** verification **/
    /*
     multiple matches:
     - will work, multiple (unexpected?) products
     one match each: {p1_fwd,  p1_rev, p1_absent} x {p2_fwd,  p2_rev, p2_absent}
     - 3 cases (fails): 1 or both absent
     - 2 cases (fails): same direction
     - 2 cases: point toward each other
     p1 = fwd, p2 = rev (fwd < rev)  or p1 = rev, p2 = fwd (fwd < rev)
     - 2 cases: point away from each other
     same as above but (fwd > rev)
     special cases: ?
     */
    /**
     * @description Determines if primers will anneal and only once (Simple PCR)
     * @param sequence
     * @param primer1
     * @param primer2
     * @returns {string|boolean} true if no error, otherwise string with error
     */
    var verifyPrimers = function verifyPrimers(sequence, primer1, primer2) {

        if (_.isEmpty(primer1) || _.isEmpty(primer2)) {
            return "a primer is not defined"
        }

        var p1 = findAnnealAllExact(sequence, primer1),
            p2 = findAnnealAllExact(sequence, primer2);

        /** check zero matches **/

        if (p1.forward.length + p1.reverse.length < 1) {
            return "primer1 : no matches"
        }
        if (p2.forward.length + p2.reverse.length < 1) {
            return "primer2 : no matches"
        }

        /** check multiple matches **/

        if (p1.forward.length + p1.reverse.length > 1) {
            return "primer1 : multiple matches"
        }
        if (p2.forward.length + p2.reverse.length > 1) {
            return "primer2 : multiple matches"
        }

        /** check directions */

        if ((p1.forward.length && p2.forward.length) || (p1.reverse.length && p2.reverse.length))
            return "primers point same direction";

        return true;

    };



    // see https://www.ncbi.nlm.nih.gov/tools/epcr/

    //wrapper function, currently only handles PCR and EIPCR
    var predict = function predict(sequence, primers) {

        //console.log(sequence, primers);

        if (primers.length != 2)
            return "Can only handle having two primers right now";

        var primer1 = primers[0],
            primer2 = primers[1];

        var verify = verifyPrimers(sequence, primer1, primer2);
        if (verify !== true)
            return verify;


        //future - not DRY
        var p1 = findAnnealAllExact(sequence, primer1),
            p2 = findAnnealAllExact(sequence, primer2);


        /** orient primers **/
        //todo - break out

        var p1pos = (!!p1.forward.length) ? +p1.forward[0] : +p1.reverse[0];
        var p2pos = (!!p2.forward.length) ? +p2.forward[0] : +p2.reverse[0];



        //pass to protocol

        //p1 forward, p2 reverse
        if (!!p1.forward.length) {
            //normal
            if (p1pos < p2pos) {
                return PCR(sequence, p1pos, (p2pos + primer2.length));
            }
            //eipcr
            else {
                return EIPCR(sequence, p1pos, (p2pos+ primer2.length));
            }
        }
        //p2 forward, p1 reverse
        else {
            //normal
            if (p1pos > p2pos) {
                return PCR(sequence, p2pos, (p1pos + primer1.length));
            }
            //eipcr
            else {
                return EIPCR(sequence, p2pos, (p1pos + primer1.length));
            }

        }
    };



    /** extension **/
    //note- currently logic for this functions exists in wrapper function

    var PCR = function PCR(sequence, forwardPrimerPos, reversePrimerPos) {
        //console.log(arguments);
        return sequence.substring(forwardPrimerPos, reversePrimerPos);
    };

    //aka Overlap Extension
    //note - currently, very simple, only two overlapping primers with defined overlap length
    var SOEing = function () {

    };

    var EIPCR = function EIPCR(sequence, forwardPrimerPos, reversePrimerPos) {
        //console.log(arguments);
        return sequence.substring(forwardPrimerPos) + sequence.substring(0, reversePrimerPos);
    };

    //30-150 bp, 2 primers w/ homology region (~20bp) and extend
    var wobble = function wobble(primer1, primer2, overlapLength) {
        //future - calculate overlapLength if not given

        return primer1 + (DNA.revcomp(primer2)).substring(overlapLength);
    };

    var PCA = function PCA() {

    };

    var RTPCR = function RTPCR() {

    };

    var Gibson = function Gibson() {
        //adds ligase (remove nicks)
        //5' exo nuclease
    };

    var SLIC = function SLIC() {
        //3' exo nuclease

    };





    /**************
     Alignment
     **************/
    //todo - fold into PCR function
    //todo - handle multiple lines
    var primerAlign = function(sequence, primers) {
        if (primers.length != 2)
            return "Can only handle having two primers right now";

        var primer1 = primers[0],
            primer2 = primers[1];

        var verify = verifyPrimers(sequence, primer1, primer2);
        if (verify !== true)
            return verify;


        //future - not DRY
        var p1 = findAnnealAllExact(sequence, primer1),
            p2 = findAnnealAllExact(sequence, primer2);


        /** orient primers **/
        //todo - break out

        var p1pos = (!!p1.forward.length) ? +p1.forward[0] : +p1.reverse[0];
        var p2pos = (!!p2.forward.length) ? +p2.forward[0] : +p2.reverse[0];

        //console.log(p1pos, p2pos, sequence, primers[0], primers[1]);

        //todo - don't assume p1 first

        var line1 = DNA.createRun(' ', p1pos);
        line1 += (!!p1.forward.length) ? primers[0] : DNA.revcomp(primers[0]);
        line1 += DNA.createRun(' ', (p2pos - p1pos - primers[0].length));
        line1 += (!!p2.forward.length) ? primers[1] : DNA.revcomp(primers[1]);
        line1 += DNA.createRun(' ', (sequence.length - line1.length));

        //note - handle line breaks outside function
        return (line1);
    };




    /**************
     ligation
     **************/



    /**
     * @description Ligates two fragments. Optionally, shows alignment of ligation
     * @param fragments {Array} NOTE Currently only two
     * @param align {boolean} Whether to show alignment (both strands)
     * @param showHTML {boolean} Include HTML tags marking original fragments and complementarity region. Recommended when displaying as HTML. Default false.
     * @param showMarks {boolean} maintain cut marks. Not relevant if showHTML. Default false.
     * @returns {string}
     */
    var ligateOld = function(fragments, align, showHTML, showMarks) {

        var fragments = parseFragments(fragments),
            blunts,
            overhangs,
            fragPair;

        if (blunts.length > 1) {
            if (blunts.length == 2) {
                console.log('2 blunt ends -- expect random products, just showing one');
                fragPair = orientFragmentsForJoin(blunts[0], blunts[1]);
            }
            else {
                //several weird products
                return 'multiple blunt end fragments -- random products'
            }
        }

        if (overhangs.length > 1) {
            if (overhangs.length == 2) {
                fragPair = orientFragmentsForJoin(overhangs[0], overhangs[1]);
                return joinTwoFragments(fragPair, align, showHTML, showMarks);
            } else {
                //multiple pairs
                console.log('more than 2 overhangs', overhangs);
            }
        }



                /* rearchitecture:

                 todo - non-terminal ends should be cut before go through

                 todo ---- each fragment should maintain both of its ends, organize pointers by fragment -> end -> pointer ----- don't separate into separate fragments


                 note - aaaaA^AGGT_Tgggg = ccccT^TCCA_Atttt


                 */



                function findConnections(overhangs) {
                    var pointers = {};

                    _.each(overhangs, function(currentHang, indout) {
                        _.each(overhangs, function(otherHang, indin) {
                            if (indout != indin) {
                                var curEnd = currentHang.end.match,
                                    otherEnd = otherHang.end.match;

                                //if ends are complimentary
                                if (curEnd == otherEnd || curEnd == DNA.revcomp(otherEnd)) {
                                    //first, check don't have the other way mapped
                                    if (!pointers[indin] || pointers[indin] != indout) {
                                        if (!pointers[indout]) {
                                            pointers[indout] = indin;
                                        } else {
                                            //points to two fragments
                                            return 'fragment ' + indout + ' has multiple complementary fragments: ' + indin + ' and ' + pointers[indout];
                                        }
                                    }
                                }
                            }
                        })
                    });

                    //can assume at this point all pointers are unique and one way
                    console.log(pointers);

                    if (!_.keys(pointers).length) return false;

                    return pointers;
                }

                function makeConnections(pointers) {
                    var products = [];
                    //todo - handle 3+ way binds (i.e. check next for partner)
                    _.each(pointers, function (value, key) {
                        fragPair = orientFragmentsForJoin(overhangs[key], overhangs[value]);
                        products.push(joinTwoFragments(fragPair, align, showHTML, showMarks))
                    });


                    console.log(products);

                    return products;
                }

                //fixme -- need to process ends again
                var inputFrags = overhangs,
                    pointers,
                    products;
                while (!!(pointers = findConnections(inputFrags)) ) {
                    products = makeConnections(pointers);
                    inputFrags = products;
                }



                return products;





                 /*var matches = {};

                 _.each(overhangs, function(overhang, outerIndex) {
                 var curHang = overhangs[outerIndex].end.match;
                 _.each(overhangs, function(overhang, InnerIndex) {
                 console.log('overhang ' + outerIndex + ' match ' + InnerIndex, curHang == overhang.end.match, curHang == DNA.revcomp(overhang.end.match), overhangs[outerIndex] !== overhang);
                 if ((curHang == overhang.end.match || curHang == DNA.revcomp(overhang.end.match)) && overhangs[outerIndex] !== overhang) {
                 matches[outerIndex] = !!matches[outerIndex] ? matches[outerIndex].push(InnerIndex) : [InnerIndex];
                 }
                 });
                 });

                 //console.log(matches);

                 //check for arrays.length > 1
                 if (_.filter(matches, function(matchArr, key) { return matchArr.length != 1; }).length)
                 return 'more than one match for some overhangs';

                 if (!matches.length) {
                 //if only two, just do the fragPair thing
                 if (_.keys(matches).length == 2) {
                 console.log('only one set of matching overhangs');

                 if (overhangs.length > _.keys(matches).length)
                 console.log('some fragments ends did not match and will be ignored');

                 //just pull one
                 var key = _.keys(matches)[0],
                 match = matches[key][0];
                 fragPair = orientFragmentsForJoin(overhangs[key], overhangs[match]);
                 }
                 else {
                 return 'multiple (' + _.keys(matches).length + ') sets of matches -- currently can only handle 2 matching overhangs';

                 }
                 } else {
                 return 'no matching overhangs'
                 }
            }
        }
        */
    };



    //DEPRECATED
    //todo - incorporate into Fragment.joinFragment
    var joinTwoFragments = function (fragPair, align, showHTML, showMarks) {

        if (!fragPair)
            return 'cut marks not defined';

        //todo - why is this here?
        if (fragPair[0].end.match != DNA.revcomp(fragPair[1].end.match))
            return 'overhangs not complimentary (5\' orientation): ' + fragPair[0].end.match + ' : ' + fragPair[1].end.match;


        var matchType = (fragPair[0].end.isBlunt) ? 'bluntmatch' : 'stickymatch';

        var finalText = '<ligate-frag>' + (fragPair[0].fragment).substring(0, fragPair[0].end.index) + '</ligate-frag>' +
            '<ligate-'+matchType+'>' +
            (showMarks ? fragPair[0].end.match : Digest.removeMarks(fragPair[0].end.match ))+
            '</ligate-'+matchType+'>' +
            (fragPair[1].fragment).substring(fragPair[1].end.index + fragPair[1].end.match.length);


        if (!!align) {
            var line2 = DNA.complement((fragPair[0].fragment).substring(0, fragPair[0].end.index)) +
                '<ligate-'+matchType+'>' +
                (showMarks ? fragPair[1].end.match : Digest.removeMarks(fragPair[1].end.match ))+
                '</ligate-'+matchType+'>' +
                '<ligate-frag>' + DNA.complement((fragPair[1].fragment).substring(fragPair[1].end.index + fragPair[1].end.match.length)) + '</ligate-frag>';

            finalText += "\n" + line2;
        }

        if (!showHTML) {
            finalText = finalText.replace(/(<([^>]+)>)/ig, '');
        }

        //fixme - handle differently -- shouldn't remove marks outside of ligation
        /*if (!showMarks && !showHTML) {
         finalText = Digest.removeMarks(finalText);
         }*/

        return finalText
    };


    //given strings, creates fragments with terminal ends
    //shouldTrim --- for internal cuts,
    //  true = sequence past overhang removed,
    //  false = split into fragments
    var parseFragments = function PCR_parseFragments (fragments, shouldTrim) {

        var parsedFrags = [];
        _.each(fragments, function(frag, indout) {

            //e.g. if pass in digest array directly, assume sorted and take first
            if (_.isArray(frag))
                frag = frag[0];

            //nonterminal marks
            if ((Digest.findOverhangs(frag, true)).length) {
                if (shouldTrim) {
                    var trimmed = Digest.trimPastInternal(frag, true);
                    parsedFrags.push(new Fragment(trimmed))

                } else {
                    _.each(Digest.makeCuts(frag), function (cutFrag) {
                        parsedFrags.push(new Fragment(cutFrag));
                    });
                }
            }
            else {
                parsedFrags.push(new Fragment(frag));
            }
        });

        return parsedFrags;
    };


    /**
     *
     * @param {Array} fragments Array of Strings WITH cut marks
     */
    //todo - add options -- HTML & alignment
    var ligate = function(fragments) {
        fragments = parseFragments(fragments, true);
        //console.log('ligate starting -- fragments:', fragments);

        //use for loop so can decrement counter
        for (var outerInd = 0; outerInd < fragments.length; outerInd++) {
            var outerFrag = fragments[outerInd];

            var toJoinIndex = outerFrag.findFirstMatchIndex(fragments);

            if (toJoinIndex >= 0) {
                //testing
                //console.log('outer frag at index ' + outerInd, outerFrag);
                //console.log('first match at index ' + toJoinIndex, fragments[toJoinIndex]);

                outerFrag.joinFragment(fragments[toJoinIndex]);

                //console.log('outer frag is now:', outerFrag, "\n\n\n\n\n");
                fragments.splice(toJoinIndex, 1);
                outerInd--;
            }
        }

        _.each(fragments, function(fragment) {
            fragment.circularize();
        });

        //console.log('LIGATE FINAL:', fragments);

        if (fragments.length == 1)
            return fragments[0].sequence;
        else
            return _.pluck(fragments, 'sequence');
    };

    return {
        //PCR
        predict : predict,

        //Annealing
        anneal : anneal,
        findAnnealSimple : findAnnealSimple,

        //Alignment
        primerAlign : primerAlign,

        //Ligate
        parseFragments : parseFragments,
        ligate : ligate


    };
}]);

Application.Dna.directive('pcrPredict', ['PCR', 'Digest', 'DNA', function(PCR, Digest, DNA) {

    return {
        restrict: 'A',
        require: 'ngModel',
        scope: {
            backbone: '=ngModel',
            primers: '='
        },
        link: function pcrPredictLink(scope, element, attr, ngModel) {
            ngModel.$render = function() {
                process();
            };

            scope.$watch(function() {
                return scope.primers[0] + scope.primers[1]
            }, function() {
                process();
            });

            function process () {
                element.text(PCR.predict(scope.backbone, scope.primers));
            }
        }
    };
}]);

Application.Dna.directive('pcrAlign', ['PCR', 'Digest', 'DNA', '$filter', function(PCR, Digest, DNA, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        scope: {
            backbone: '=ngModel',
            primers: '='
        },
        link: function pcrPredictLink(scope, element, attr, ngModel) {
            ngModel.$render = function() {
                process();
            };

            scope.$watch('primers', process, true);

            function process () {
                var alignment = PCR.primerAlign(scope.backbone, scope.primers);

                //console.log(alignment);

                //todo -- pass in width of element to breaklines
                var charNum = 57;

                alignment = $filter('breakLines')(alignment, charNum, "*").split('*');

                var backboneText = $filter('breakLines')(scope.backbone, charNum, "*").split('*');

                //console.log(alignment, backboneText);

                var finalText = "";
                for (var i = 0; i < backboneText.length; i++) {
                    finalText += alignment[i] + "\n" + backboneText[i] + "\n";
                }

                //console.log(finalText);
                element.html(finalText)
            }



        }
    }
}]);

Application.Dna.directive('ligateAlign', ['PCR', 'Digest', 'DNA', '$compile', '$filter', function(PCR, Digest, DNA, $compile, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        scope: {
            fragments: '=ngModel'
        },
        link: function pcrPredictLink(scope, element, attr, ngModel) {
            ngModel.$render = function() {
                process();
            };

            scope.$watch('fragments', function() {
                process();
            }, true);

            function process () {
                var alignment = PCR.ligate(scope.fragments, true, true);
                console.log(alignment);

                if (_.isArray(alignment)) {
                    alignment = "did not ligate to completion : " + JSON.stringify(alignment);
                } else {
                    alignment = $filter('DigestCuts')(alignment, true);
                }

                element.html($compile('<span>' + alignment + '</span>')(scope))
            }
        }
    }
}]);

Application.Dna.directive('ligateFrag', function() {
    return {
        restrict: 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="Initial Fragment" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        link: function(scope, element, attrs) {
            element.css('color', '#faa');
        }
    }
});

Application.Dna.directive('ligateStickymatch', function() {
    return {
        restrict: 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="Sticky-end Complementary Region" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        link: function(scope, element, attrs) {
            element.css('color', '#6b6');
        }
    }
});

Application.Dna.directive('ligateBluntmatch', function() {
    return {
        restrict: 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="Note! Blunt ends will only yield this product 50% of the time (fragment direction is random)" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        link: function(scope, element, attrs) {
            element.css('color', '#6b6');
        }
    }
});
