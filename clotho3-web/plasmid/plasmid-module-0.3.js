'use strict';

/*

 todo
 - restrict input to one of multiple sets
    - ACGT
    - plus others (degenerates - e.g. N, Y, R)
    - RNA
 - revcomp method
 - other parameters
    - circular
    - linear
    - GC content / emlting point (from server)
 - start / stop codons
 - check both directions
 - add feature directive -- tooltip?

 how will funciton
 - use uniquely named tags with feature directive attr (function as tooltip)
 - features - rollover to only show certain feature (easy - just strip all other tags)

 algorithm
 INITIAL
 - create locations hash mapping pos -> tag
 - start from back, add in tags
 - future - check reverse direction too

 CHANGES
 - add feature by highlight
 - just add to HTML
 - add sequence
 - just add to sequence
 - add regexp / feature indirectly
 - requires recompile
 - pull out tags and locations from HTML (start at front)
 - search for regexp in text, add locations
 */


Application.Plasmid.controller('PlasmidCtrl', ['$scope', '$window', '$document', function($scope, $window, $document) {
    $scope.sequence = "GATCTgttgacggctaGCTCAGTCCTAGGTagctACAGTGCTAGCTCTCTGGAGATTAACGAGGAGAAATACTAGATGGTTCATGATCATAAgcttgaattagccaaacttattcgcaactatgagacgaatagaaaagaatgtctaaattccagatataatgaaacacttttacgaagtgattatcttgatccattttttgaacttcttggctgggatattaaaaataaagctggaaaaccgactaatgaaagagaggttgtcttggaagaggcacttaaagcaagtgcatctgaacattctaaaaaaccagattatacattcagacttttttctgaaagaaagtttttcttggaagctaaaaaaccatcagttcatattgaatcggataatgaaactgctaaacaagtgcgaagatatggctttaccgccaaactaaaaatttcagttttatcaaattttgaatatttagttatttatgatacctctgtaaaggttgatggtgatgatacctttaataaggcacgtataaaaaaataccattacacagagtatgaaactcactttgatgaaatttgtgacttattaggaagagagtccgtttactctgggaattttgataaagaatggttgagtatcgaaaataaaattaatcacttttctgtagataccttatttttaaaacagattaatacatggcgtctattgcttggtgaagaaatctataagtatcaacctacgatacaagagaatgagcttaatgacattgtacagagctatctgaatagaattAGATCTatttttttgagagtctgtgaagatagaaatttagagacttatcagacattactgaattttgcttcaagtaatgatttctccgctcttattgataagtttaagcaggcagatcgttgctataattcaggcctatttgatcaattgcttacagagcaaattattgaggatattagttctgtattttgggtaatcattaagcaattatattatccagaaagtccttattcatttagtgtgttctcttcggatattttaggtaatatttacgaaatatttttatctgagaaattagtaattaatcaaagcagagttgagttagtcaagaaaccagagaatttagatagagacattgtcacaacaccaacctttattattaatgacatcttgagaaatacggttctaccgaagtgctatggaaaaacagatatagaaattctacagctaaaatttgctgatattgcttgtggttcgggagcatttttactggagttgttccaattacttaatgatactctagttgactattatttaagtagtgatacttctcaattaattccaacaggtatcggtacttataagctgtcttatgaaatcaagagaaaggttctattaagttgtatttttggcatagataaggacttaaatgctgtagaggctgcaaagttcggattgttgctaaaattattagagggtgaagacgtacaatctatagctaatattagaccagttctcccagatttattagataacatactttttggtaacagtttattagaaccagaaaaagtcgagcttgatcatcaggtagaagtaaatccgttagatttttcGGATCTGaaAtttgatgtaattgttggcaaccctccatatatgaaatcagaggatatgaagaatattactcctttggagttacctttatataagaaaaactatgtttctgcttataagcaatttgataaatatttcttgttcttagagcggggtttagctctattaaaagaagagggaatacttggatatattgttccaagtaaatttactaaagtgggtgcagggaaaaagttacgggaattactaacagataagggttatcttgactctattgtttcttttggtgctaatcaaatatttcaggataaaacaacttatacttgtttacttattttaagaaaaactccAcatactgattttaaatatgcagaggttcgtaatttaattgactggaaagtgcgtaaagctgatgctatggaattttcctctcaacaactgagtacattgcaaagtgatgcgtggattttaattccatctgaattaatctcagtttatcatcagatattagcacaaagccaaaagctagaggatattgtcggtattgataatatatttaatgggattcaaaccagtgctaatgatgtctatatttttgtgccaactcatgaggatactgaaaactattattttataaagaaaggacaagagtacaaaattgaaaaggaaattacgaagccttattttaaaacaacgagtggtgaggataacttatatacttaccgtactttcaagcctaatgcccgagtcatttatccgtatactcaaactgagagtagtgtagaactaattcctttagatgaaatacgagaaatttttcctttagcatacaaatatttaatgtcgcttaagttcgttttaagtagccccaaacgagatataaaacctagacctaaaacaacaaatgaatggcataggtatggacggcatcaaagtctcgataattgtgggttgagtcagaaaattattgtaggtgtgctttcagttggtgataagtacgctatagatacttatggaacgttgatttcatcaggcggtacggctggatactgtgtggttgctcttccagatgattgtaaatattcaatttattatttacaggcaattttaaactcaaaatatttagagtggtttagtgccttacatggagaagttttccgaggtggttatattgctaggggaactaaggtgcttaagaacttgcctattaggaaaattgattttgataatcttgaagaagcaaatctacatgatctaattgcgaccaagcaaaaagagcttatagagatttatgacaaaatagatgttaatgtaaataataaaagagttctgaccccattgcaacgtatgtttaaacgagagaaagaggttttagaccaattgttgagtcgactgtataacttaggtgtagatgattccttgatcccttatattaaggatttgtatgaagctcattaaGGATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAATAATTGCTCATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTGCCATCCAGCTGATATCCCCTATAGTGAGTCGTATTACATGGTCATAGCTGTTTCCTGGCAGCTCTGGCCCGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGCCTCCTCTAGACCAGCCAGGACAGAAATGCCTCGACTTCGCTGCTGCCCAAGGTTGCCGGGTGACGCACACCGTGGAAACGGATGAAGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGGGGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAAACATCATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATCGAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGCGGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAAACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAACATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAGGATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCTGGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGCAAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTATCAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCCTCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTAGTCGGCAAATAACCCTCGAGCCACCcatgaccaaaatcccttaacgGCATGCgcaccgccggacatcagcgctagcggagtgtatactggcttactatgttggcactgatgagggtgtcagtgaagtgcttcatgtggcaggagaaaaaaggctgcaccggtgcgtcagcagaatatgtgatacaggatatattccgcttcctcgctcactgactcgctacgctcggtcgttcgactgcggcgagcggaaatggcttacgaacggggcggagatttcctggaagatgccaggaagatacttaacagggaagtgagagggccgcggcaaagccgtttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaagaagatcatcttattaatcagataaaatatttctagAGGCCTcccctgattctgtggataaccGTcctaggTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctggaattcatgAtacgaa";

    $scope.editable = true;
    $scope.textSelected = false;

    $scope.featureList = [
        {
            "label" : "pCon Promoter",
            "pos" : {
                "start" : 40,
                end : 46
            },
            "css" : {
                "color" : "",
                "background" : "#99EEEE"
            }
        },
        {
            "label" : "Another Promoter",
            "pos" : {
                "start" : 90,
                end : 106
            },
            "css" : {
                "color" : "",
                "background" : "#EE99EE"
            }
        },
        {
            "label" : "Binding Site Promoter",
            "pos" : {
                "start" : 120,
                end : 146
            },
            "css" : {
                "color" : "",
                "background" : "#EEEE99"
            }
        },
        {
            "label" : "Xbal II",
            "match" : "TAGCTA",
            "css" : {
                "color" : "#FF88CC",
                "background" : ""
            }
        },
        {
            "label" : "Feature1",
            "match" : "ACAGTG",
            "css" : {
                "color" : "#BBDD77",
                "background" : ""
            }
        },
        {
            "label" : "Feature2",
            "match" : "ATCGAT",
            "css" : {
                "color" : "#33DDDD",
                "background" : ""
            }
        },
        {
            "label" : "Feature3",
            "match" : "attaagcaa",
            "css" : {
                "color" : "#DD44DD",
                "background" : ""
            }
        },
        {
            "label" : "Feature4",
            "match" : "TAATTTCAT",
            "css" : {
                "color" : "#FF8811",
                "background" : ""
            }
        },
        {
            "label" : "Feature5",
            "match" : "AGATCT",
            "css" : {
                "color" : "#CC22BB",
                "background" : ""
            }
        },
        {
            "label" : "Feature6",
            "match" : "TACGAA",
            "css" : {
                "color" : "#22FF88",
                "background" : ""
            }
        },
        {
            "label" : "Feature7",
            "match" : "TATGGAA",
            "css" : {
                "color" : "#FFBB33",
                "background" : ""
            }
        },
        {
            "label" : "Restriction Site",
            "match" : "TAGCTAG",
            "css" : {
                "color" : "",
                "background" : "#EECC99"
            }
        }
    ];



    //todo - expose, make editable
    $scope.reg_match = /^[acgtACGT]+$/;
    $scope.reg_color = /^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/;
    $scope.reg_pos = /^[0-9]+\-[0-9]+$/;


    $scope.emptyFeat = function() {
        return {label: "", pos: '', match: '', css : {color: '', background: '' }};
    };
    
    $scope.logFeatures = function() {console.log($scope.featureList)};

    $scope.featureValid = function(feat) {
        return (
                (feat.match != '' && angular.isDefined(feat.match)) ?
                $scope.reg_match.test(feat.match) :
                (angular.isDefined(feat.pos) && $scope.reg_pos.test(feat.pos))
            ) &&
            ( (angular.isDefined(feat.css.color) && feat.css.color != '') ?
                ($scope.reg_color.test(feat.css.color)) :
                (angular.isDefined(feat.css.background) && ($scope.reg_color.test(feat.css.background))))
    };

    $scope.addFeature = function(feat) {
        if ($scope.featureValid(feat)) {

            if (!feat.label) {
                feat.label = "New Feature" + Date.now();
            }

            if (feat.match)
                delete feat['pos'];
            else {
                var positions = feat.pos.match(/^([0-9]+)\-([0-9]+)$/);
                feat.pos.start = positions[1];
                feat.pos.end = positions[2];
            }

            if (feat.css.color == '' && feat.css.background == '')
                feat.css.background = '#'+Math.floor(Math.random()*16777215).toString(16);

            console.log(feat);

            $scope.featureList.push(feat);
            $scope.new = $scope.emptyFeat();
        } else {
            console.log("invalid");
        }
    };


    //todo - checks:
    // select opp direction
    // pos is for whole string (across text nodes)
    $scope.addFeatureSelection = function () {

        var feature = $scope.emptyFeat();
        if ($window.getSelection) {
            var sel = $window.getSelection();
            if (typeof sel != 'undefined' && sel.rangeCount) {
                var container = angular.element("<div>");
                for (var i = 0, len = sel.rangeCount; i < len; ++i) {
                    container.append(sel.getRangeAt(i).cloneContents());
                }
                console.log(sel.getRangeAt(0).cloneRange());
                //feature = container[0].innerHTML;

                console.log(sel);
                feature.pos = {"start" : sel.baseOffset, "end" : sel.extentOffset}


                feature.css.background = '#'+Math.floor(Math.random()*16777215).toString(16);

            }
        }
        else if (typeof $document.selection != "undefined") { //IE -- todo doesn't really work
            if ($document.selection.type == "Text") {
                feature = $document.selection.createRange().htmlText;
            }
        }
        console.log(feature);
    };

}]);

Application.Plasmid.directive('plasmidEditor', ['$parse', '$timeout', '$filter', '$compile', '$document', function($parse, $timeout, $filter, $compile, $document) {

    return {
        restrict: "A",
        require:'ngModel',
        //rep1ace: true,
        scope : {
            features: '=',
            editable: '@',
            textSelected: '=', // todo - move into this directive
            sequence: '=ngModel'
        },
        link: function(scope, element, attrs, ngModel) {
            attrs.$set('spellcheck', false);
            attrs.$set('contenteditable', true); //todo - check controller

            //scope.features = $parse(attrs.features)(scope) || [];
            //scope.textSelected = $parse(attrs.textSelected)(scope);

            console.log(ngModel);

            /* key functions  */

            function genFiltered (text) {
                text = typeof text != 'undefined' ? text : ngModel.$viewValue;
                return $filter('features')(text, scope.features);
            }

            /* updating view */

            function setOutput (html) {
                //var x = $compile(html)(scope);
                //console.log(x);
                element[0].innerHTML = html;
                $compile(element.contents())(scope);
            }


            //model -> view

            /*scope.$watch(function () {
                return ngModel.$modelValue;
            }, function (modelValue) {
                console.log('model change');
                //setOutput(genFiltered())
            });*/


            //fixme - why isn't render being called with updates??
            ngModel.$render = function() {
                //console.log(ngModel.$viewValue);
                console.log('render');
                element.html(genFiltered(ngModel.$viewValue));
            };

            scope.$watch(function() {
                return ngModel.$viewValue;
            }, function(viewVal) {
                console.log('view watch');
                //element.html(genFiltered(viewVal));
            });

            scope.$watch(function() {
                return ngModel.$modelValue;
            }, function(modelVal) {
                console.log('model watch');
                //element.html(genFiltered(viewVal));
            });


            /* watchers */
            scope.$watchCollection('features', function(newVal, oldVal) {
                if (!!newVal && !!oldVal) {
                    setOutput(genFiltered());
                }
            });

            element.bind( "keyup blur change" , function() {
                console.log('change to view model');
                scope.$apply(
                    ngModel.$setViewValue(element.html())
                );
            });

            element.bind('mouseup', function() {
                if (typeof window.getSelection != "undefined") {
                    scope.textSelected = !!window.getSelection().toString();
                } else if (typeof document.selection != "undefined" && document.selection.type == "Text") {
                    scope.textSelected = !!document.selection.createRange().text;
                }
            });
        }
    }
}]);

Application.Plasmid.filter('features', [function() {
    return function (text, features) {
        if (features && angular.isArray(features) && angular.isString(text)) {

            console.log(features);

            var html = text;
            text = html.replace(/(<([^>]+)>)/ig, "");
            
            //console.log(html);
            //console.log(text);


            //future - pull tags out of HTML and save locations

            //var str=html.replace(/ /ig,'<span>$3</span>');
            //console.log(str + '\n\n\n\n\n\n\n\n\n\n\n\n');


            //create location map
            var locations = {};
            angular.forEach(features, function(feat, featIndex) {
                if (feat.match) {
                    for (var index, offset = 0, search = angular.lowercase(text);
                         (index = search.indexOf(angular.lowercase(feat.match), offset)) > -1;
                         offset = index + feat.match.length
                        ) {
                        //start
                        var start = index;
                        locations[start] ?
                            locations[start]['start'].push(featIndex) :
                            locations[start] = {"start" : [featIndex], "end" : []};
                        //end
                        var end = start + feat.match.length;
                        locations[end] ?
                            locations[end]['end'].push(featIndex) :
                            locations[end] = {"start" : [], "end" : [featIndex]};

                    }
                } else {
                    //start
                    locations[feat.pos.start] ?
                        locations[feat.pos.start]['start'].push(featIndex) :
                        locations[feat.pos.start] = {"start" : [featIndex], "end" : []};
                    //end
                    locations[feat.pos.end] ?
                        locations[feat.pos.end]['end'].push(featIndex) :
                        locations[feat.pos.end] = {"start" : [], "end" : [featIndex]};
                }
            });
            //console.log(locations);



            //future - check reverse direction too



            //loop from end, saving tag locations
            var reversed = [];
            angular.forEach(locations, function(key, val) {
                reversed.unshift(val);
            });



            var newText = text,
                backlog = [];

            /* todo need to deal with scenarios like this:
            //note - angular will add class="ng-scope" to tag
             <1> xx    <2> xx         <3> xx     </2> xx          </3> xx      </1>      ->
             <1> xx</1><1-2> xx </1-2><1-2-3> xx </1-2-3><1-3> xx </1-3><1> xx </1>
             */


            angular.forEach(reversed, function(index) {
                angular.forEach(locations[index], function(value, key) {
                    if (text.length < index) {
                        console.log("string too short");
                        return;
                    }

                    angular.forEach(value, function(featIndex) {
                        var feat = features[featIndex];
                        var featName = angular.lowercase(feat.label).replace(/[ _]/gi, '');

                        //console.log(feat);
                        //console.log(key + " " + feat.label);

                        if (key == 'start') {
                            newText = newText.slice(0, index) +
                                '<feat' + featName + ' ' +
                                'feature="' + feat.label + '" ' +
                                'style="' +
                                (feat.css.color ? "color: " + feat.css.color + ";" : "") +
                                (feat.css.background ? "background-color: " + feat.css.background : "") +
                                '">' +
                                newText.slice(index);
                        } else {
                            newText = newText.slice(0, index) +
                                '</feat' + featName + '>' +
                                newText.slice(index);
                        }
                    });

                })
            });
            //note - text is filled with features at appropriate locations. tags might overlap.
            //console.log(newText);













            /**** REGEXPS *****/
            var findAllTags = /<\/?(feat[A-Z0-9]*)\b[^>]*>/gi;

            var findFeatsInclusive = /<(feat[A-Z0-9]*)\b[^>]*>.*?<\/\1>/gi;

            var findSingleOverlaps = /<feat([A-Z0-9]*)\b[^>]*>.*?<feat([A-Z0-9]*)\b[^>]*>[^\/]*?(?!(\/feat\2))\/feat\1>/gi;

            //note - find features with tag inside (overlapping)
            //todo - fix
            var findOverlaps = /<(feat[A-Z0-9]*)\b[^>]*>.*?<(feat[A-Z0-9]*)\b[^>]*>.*?(?!\2).*?<\/\1>/gi;

            //todo - fix
            var findSingleNested = /<(feat[A-Z0-9]*)\b[^>]*>.*?<(feat[A-Z0-9]*)\b[^>]*>.*?<\/\2>.*?<\/\1>/gi;



            var feature_reg = /<feat([A-Z0-9]*) feature="(.*?)" ([^>]*)>(.*?)<feat([A-Z0-9]*) feature="(.*?)" ([^>]*)>([^\/]*?)<(?!(\/feat\3))\/feat\1>/gi;
            var feature_replacer = function(match, f1, n1, c1, s1, f2, n2, c2, s2, ignore, index){
                /*console.log(arguments);
                console.log(index);
                console.log(match);*/

                var string = '<feat' + f1 + ' feature="' + n1 + '" ' + c1 + '>' + s1 +
                    '</feat' + f1 + '>' +
                    '<feat' + f1 + '-' +  f2 + ' feature="' + n1 + ', ' + n2 + '" ' + c1 + '>' + s2 +
                    '</feat' + f1 + '-' +  f2 + '>' +
                    '<feat' + f2 + ' feature="' + n2 + '" ' + c2 + '>';

                return string;

            };

            var overlap;
            while (overlap = findSingleOverlaps.exec(newText)) {
                //var inner = reg.exec(overlap[0]);
                //console.log(inner);

                var corrected = overlap[0].replace(feature_reg, feature_replacer);

                newText = newText.replace(overlap[0], corrected);
            }

            //console.log(newText);







            return newText;
        } else {
            return text;
        }
    };
}]);

//todo - interaction with Plasmid Service? or just tooltip?
Application.Plasmid.directive('feature', ['$tooltip', function($tooltip) {

    return {
        restrict : 'A',
        replace: true,
        scope: {
            feature: '@'
        },
        transclude:true,
        template: '<feature><span tooltip="{{ feature }}" append-to-body ng-transclude></span></feature>', //todo - get append-to-body working so font works
        link: function(scope, element, attrs, ctrl) {
            //borrow from angularUI tooltip?
            //attrs.$set('style', "background-color: #FF0000");
            
//          console.log(scope);

            element.bind('mouseenter', function() {
                console.log(scope.feature);
            });

            //todo - simple notification - ask if want to split or what
           /*
           scope.$watch(function() {
                return element.text();
            }, function( newval, oldval) {
                if (!!newval && !!oldval && newval != oldval) {
                    alert('changing a feature!')
                }
            });
            */

        }
    }
}]);