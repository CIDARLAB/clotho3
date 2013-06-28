'use strict';

/*

 todo
 - restrict input to one of two sets
 - ACGT
 - plus others (e.g. N, Y, R)
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


Application.Plasmid.controller('PlasmidCtrl', ['$scope', function($scope) {
    $scope.sequence = "GATCTgttgacggctaGCTCAGTCCTAGGTagctACAGTGCTAGCTCTCTGGAGATTAACGAGGAGAAATACTAGATGGTTCATGATCATAAgcttgaattagccaaacttattcgcaactatgagacgaatagaaaagaatgtctaaattccagatataatgaaacacttttacgaagtgattatcttgatccattttttgaacttcttggctgggatattaaaaataaagctggaaaaccgactaatgaaagagaggttgtcttggaagaggcacttaaagcaagtgcatctgaacattctaaaaaaccagattatacattcagacttttttctgaaagaaagtttttcttggaagctaaaaaaccatcagttcatattgaatcggataatgaaactgctaaacaagtgcgaagatatggctttaccgccaaactaaaaatttcagttttatcaaattttgaatatttagttatttatgatacctctgtaaaggttgatggtgatgatacctttaataaggcacgtataaaaaaataccattacacagagtatgaaactcactttgatgaaatttgtgacttattaggaagagagtccgtttactctgggaattttgataaagaatggttgagtatcgaaaataaaattaatcacttttctgtagataccttatttttaaaacagattaatacatggcgtctattgcttggtgaagaaatctataagtatcaacctacgatacaagagaatgagcttaatgacattgtacagagctatctgaatagaattAGATCTatttttttgagagtctgtgaagatagaaatttagagacttatcagacattactgaattttgcttcaagtaatgatttctccgctcttattgataagtttaagcaggcagatcgttgctataattcaggcctatttgatcaattgcttacagagcaaattattgaggatattagttctgtattttgggtaatcattaagcaattatattatccagaaagtccttattcatttagtgtgttctcttcggatattttaggtaatatttacgaaatatttttatctgagaaattagtaattaatcaaagcagagttgagttagtcaagaaaccagagaatttagatagagacattgtcacaacaccaacctttattattaatgacatcttgagaaatacggttctaccgaagtgctatggaaaaacagatatagaaattctacagctaaaatttgctgatattgcttgtggttcgggagcatttttactggagttgttccaattacttaatgatactctagttgactattatttaagtagtgatacttctcaattaattccaacaggtatcggtacttataagctgtcttatgaaatcaagagaaaggttctattaagttgtatttttggcatagataaggacttaaatgctgtagaggctgcaaagttcggattgttgctaaaattattagagggtgaagacgtacaatctatagctaatattagaccagttctcccagatttattagataacatactttttggtaacagtttattagaaccagaaaaagtcgagcttgatcatcaggtagaagtaaatccgttagatttttcGGATCTGaaAtttgatgtaattgttggcaaccctccatatatgaaatcagaggatatgaagaatattactcctttggagttacctttatataagaaaaactatgtttctgcttataagcaatttgataaatatttcttgttcttagagcggggtttagctctattaaaagaagagggaatacttggatatattgttccaagtaaatttactaaagtgggtgcagggaaaaagttacgggaattactaacagataagggttatcttgactctattgtttcttttggtgctaatcaaatatttcaggataaaacaacttatacttgtttacttattttaagaaaaactccAcatactgattttaaatatgcagaggttcgtaatttaattgactggaaagtgcgtaaagctgatgctatggaattttcctctcaacaactgagtacattgcaaagtgatgcgtggattttaattccatctgaattaatctcagtttatcatcagatattagcacaaagccaaaagctagaggatattgtcggtattgataatatatttaatgggattcaaaccagtgctaatgatgtctatatttttgtgccaactcatgaggatactgaaaactattattttataaagaaaggacaagagtacaaaattgaaaaggaaattacgaagccttattttaaaacaacgagtggtgaggataacttatatacttaccgtactttcaagcctaatgcccgagtcatttatccgtatactcaaactgagagtagtgtagaactaattcctttagatgaaatacgagaaatttttcctttagcatacaaatatttaatgtcgcttaagttcgttttaagtagccccaaacgagatataaaacctagacctaaaacaacaaatgaatggcataggtatggacggcatcaaagtctcgataattgtgggttgagtcagaaaattattgtaggtgtgctttcagttggtgataagtacgctatagatacttatggaacgttgatttcatcaggcggtacggctggatactgtgtggttgctcttccagatgattgtaaatattcaatttattatttacaggcaattttaaactcaaaatatttagagtggtttagtgccttacatggagaagttttccgaggtggttatattgctaggggaactaaggtgcttaagaacttgcctattaggaaaattgattttgataatcttgaagaagcaaatctacatgatctaattgcgaccaagcaaaaagagcttatagagatttatgacaaaatagatgttaatgtaaataataaaagagttctgaccccattgcaacgtatgtttaaacgagagaaagaggttttagaccaattgttgagtcgactgtataacttaggtgtagatgattccttgatcccttatattaaggatttgtatgaagctcattaaGGATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAATAATTGCTCATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTGCCATCCAGCTGATATCCCCTATAGTGAGTCGTATTACATGGTCATAGCTGTTTCCTGGCAGCTCTGGCCCGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGCCTCCTCTAGACCAGCCAGGACAGAAATGCCTCGACTTCGCTGCTGCCCAAGGTTGCCGGGTGACGCACACCGTGGAAACGGATGAAGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGGGGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAAACATCATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATCGAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGCGGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAAACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAACATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAGGATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCTGGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGCAAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTATCAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCCTCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTAGTCGGCAAATAACCCTCGAGCCACCcatgaccaaaatcccttaacgGCATGCgcaccgccggacatcagcgctagcggagtgtatactggcttactatgttggcactgatgagggtgtcagtgaagtgcttcatgtggcaggagaaaaaaggctgcaccggtgcgtcagcagaatatgtgatacaggatatattccgcttcctcgctcactgactcgctacgctcggtcgttcgactgcggcgagcggaaatggcttacgaacggggcggagatttcctggaagatgccaggaagatacttaacagggaagtgagagggccgcggcaaagccgtttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaagaagatcatcttattaatcagataaaatatttctagAGGCCTcccctgattctgtggataaccGTcctaggTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctggaattcatgAtacgaa";

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
            console.log(feat);

            $scope.featureList.push(feat);
            $scope.new = $scope.emptyFeat();
        } else {
            console.log("invalid");
        }
    }
}]);

Application.Plasmid.directive('plasmidEditor', ['$parse', '$timeout', '$filter', '$compile', '$document', function($parse, $timeout, $filter, $compile, $document) {

    return {
        restrict: "A",
        require:'ngModel',
        //rep1ace: true,
        link: function(scope, element, attrs, ngModel) {
            attrs.$set('spellcheck', false);
            attrs.$set('contenteditable', true);

            scope.features = $parse(attrs.features)(scope) || [];

            /*

            //generate the mirror

            var parent = element.parent(),
                output = angular.element('<pre>'),
                label = angular.element('<label>');

            parent.append(output);
            parent[0].className = 'ldt ' + element[0].className;


            //OLD VERSION
            //parent.append(label);
            //element[0].parentNode.replaceChild(parent, element);
            //label.append(element);


            element.bind('mouseover', function(e) {
                element.css('pointer-events', 'none');
            });
            output.bind('click', function(e) {
                e.preventDefault();
                element.css('pointer-events', 'auto')
            });

            */

            /* key functions  */

            function genFiltered (text) {
                text = typeof text != 'undefined' ? text : ngModel.$modelValue;
                return $filter('features')(text, scope.features);
            }

            /* updating view */

            function setOutput (html) {
                //var x = $compile(html)(scope);
                //console.log(x);
                element[0].innerHTML = html;
                $compile(element.contents())(scope);
            }

            function setViewValue(html) {
                ngModel.$setViewValue(html);
            }

            /* watchers */
            scope.$watchCollection('features', function(newVal, oldVal) {
                if (!!newVal && !!oldVal) {
                    setOutput(genFiltered());
                }
            });

            //model -> view
            scope.$watch(function () {
                return ngModel.$modelValue;
            }, function (modelValue) {
                setOutput(genFiltered())
            });

            //view -> model
            // detect all changes to the textarea,
            // including keyboard input, cut/copy/paste, drag & drop, etc
            if( element.addEventListener ){
                // standards browsers: oninput event
                element.bind( "input", function() {
                    scope.$apply(function () {
                        setViewValue(genFiltered());
                    });
                }, false );
            } else {
                // MSIE: detect changes to the 'value' property
                element.bind( "onpropertychange",
                    function(e){
                        if( e.propertyName.toLowerCase() === 'value' ){
                            scope.$apply(function () {
                                setViewValue(genFiltered());
                            });
                        }
                    }
                );
            }
        }
    }
}]);

Application.Plasmid.filter('features', [function() {
    return function (text, features) {
        if (features && angular.isArray(features) && angular.isString(text)) {
            var html = text;
            text = html.replace(/(<([^>]+)>)/ig, html);
            
            console.log(html);
            console.log(text);


            //future - pull tags out of HTML and save locations

            //var str=html.replace(/ /ig,'<span>$3</span>');
            //console.log(str + '\n\n\n\n\n\n\n\n\n\n\n\n');


            //create location map
            var locations = {};
            angular.forEach(features, function(feat, ing) {
                if (feat.match) {
                    for (var index, offset = 0, search = angular.lowercase(text);
                         (index = search.indexOf(angular.lowercase(feat.match), offset)) > -1;
                         offset = index + feat.match.length
                        ) {
                        //start
                        var start = index;
                        locations[start] ?
                            locations[start]['start'].push(feat) :
                            locations[start] = {"start" : [feat], "end" : []};
                        //end
                        var end = start + feat.match.length;
                        locations[end] ?
                            locations[end]['end'].push(feat) :
                            locations[end] = {"start" : [], "end" : [feat]};

                    }
                } else {
                    //start
                    locations[feat.pos.start] ?
                        locations[feat.pos.start]['start'].push(feat) :
                        locations[feat.pos.start] = {"start" : [feat], "end" : []};
                    //end
                    locations[feat.pos.end] ?
                        locations[feat.pos.end]['end'].push(feat) :
                        locations[feat.pos.end] = {"start" : [], "end" : [feat]};
                }
            });
            //console.log(locations);



            //future - check reverse direction too



            //loop from end, inserting tags
            var reversed = [];
            angular.forEach(locations, function(key, val) {
                reversed.unshift(val);
            });
            //console.log(reversed);


            var newText = text;
            angular.forEach(reversed, function(index) {
                angular.forEach(locations[index], function(value, key) {
                    if (text.length < index) {
                        console.log("string too short");
                        return;
                    }

                    //todo - check last tag, see if closing

                    angular.forEach(value, function(feat) {
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

            console.log(newText);


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

            console.log(newText);







            return newText;
        } else {
            return text;
        }
    };
}]);

//todo - interaction with Plasmid Service? or just tooltip?
Application.Plasmid.directive('feature', [function() {
    return {
        restrict : 'A',
        link: function(scope, element, attrs, controller) {
            //borrow from angularUI tooltip?
            //attrs.$set('style', "background-color: #FF0000");

            element.on('mouseenter', function() {
                console.log(attrs.feature);
            })

        }
    }
}]);