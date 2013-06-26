'use strict';

Application.Plasmid.controller('PlasmidCtrl', ['$scope', function($scope) {
    $scope.sequence = "GATCTgttgacggctaGCTCAGTCCTAGGTACAGTGCTAGCTCTCTGGAGATTAACGAGGAGAAATACTAGATGGTTCATGATCATAAgcttgaattagccaaacttattcgcaactatgagacgaatagaaaagaatgtctaaattccagatataatgaaacacttttacgaagtgattatcttgatccattttttgaacttcttggctgggatattaaaaataaagctggaaaaccgactaatgaaagagaggttgtcttggaagaggcacttaaagcaagtgcatctgaacattctaaaaaaccagattatacattcagacttttttctgaaagaaagtttttcttggaagctaaaaaaccatcagttcatattgaatcggataatgaaactgctaaacaagtgcgaagatatggctttaccgccaaactaaaaatttcagttttatcaaattttgaatatttagttatttatgatacctctgtaaaggttgatggtgatgatacctttaataaggcacgtataaaaaaataccattacacagagtatgaaactcactttgatgaaatttgtgacttattaggaagagagtccgtttactctgggaattttgataaagaatggttgagtatcgaaaataaaattaatcacttttctgtagataccttatttttaaaacagattaatacatggcgtctattgcttggtgaagaaatctataagtatcaacctacgatacaagagaatgagcttaatgacattgtacagagctatctgaatagaattAGATCTatttttttgagagtctgtgaagatagaaatttagagacttatcagacattactgaattttgcttcaagtaatgatttctccgctcttattgataagtttaagcaggcagatcgttgctataattcaggcctatttgatcaattgcttacagagcaaattattgaggatattagttctgtattttgggtaatcattaagcaattatattatccagaaagtccttattcatttagtgtgttctcttcggatattttaggtaatatttacgaaatatttttatctgagaaattagtaattaatcaaagcagagttgagttagtcaagaaaccagagaatttagatagagacattgtcacaacaccaacctttattattaatgacatcttgagaaatacggttctaccgaagtgctatggaaaaacagatatagaaattctacagctaaaatttgctgatattgcttgtggttcgggagcatttttactggagttgttccaattacttaatgatactctagttgactattatttaagtagtgatacttctcaattaattccaacaggtatcggtacttataagctgtcttatgaaatcaagagaaaggttctattaagttgtatttttggcatagataaggacttaaatgctgtagaggctgcaaagttcggattgttgctaaaattattagagggtgaagacgtacaatctatagctaatattagaccagttctcccagatttattagataacatactttttggtaacagtttattagaaccagaaaaagtcgagcttgatcatcaggtagaagtaaatccgttagatttttcGGATCTGaaAtttgatgtaattgttggcaaccctccatatatgaaatcagaggatatgaagaatattactcctttggagttacctttatataagaaaaactatgtttctgcttataagcaatttgataaatatttcttgttcttagagcggggtttagctctattaaaagaagagggaatacttggatatattgttccaagtaaatttactaaagtgggtgcagggaaaaagttacgggaattactaacagataagggttatcttgactctattgtttcttttggtgctaatcaaatatttcaggataaaacaacttatacttgtttacttattttaagaaaaactccAcatactgattttaaatatgcagaggttcgtaatttaattgactggaaagtgcgtaaagctgatgctatggaattttcctctcaacaactgagtacattgcaaagtgatgcgtggattttaattccatctgaattaatctcagtttatcatcagatattagcacaaagccaaaagctagaggatattgtcggtattgataatatatttaatgggattcaaaccagtgctaatgatgtctatatttttgtgccaactcatgaggatactgaaaactattattttataaagaaaggacaagagtacaaaattgaaaaggaaattacgaagccttattttaaaacaacgagtggtgaggataacttatatacttaccgtactttcaagcctaatgcccgagtcatttatccgtatactcaaactgagagtagtgtagaactaattcctttagatgaaatacgagaaatttttcctttagcatacaaatatttaatgtcgcttaagttcgttttaagtagccccaaacgagatataaaacctagacctaaaacaacaaatgaatggcataggtatggacggcatcaaagtctcgataattgtgggttgagtcagaaaattattgtaggtgtgctttcagttggtgataagtacgctatagatacttatggaacgttgatttcatcaggcggtacggctggatactgtgtggttgctcttccagatgattgtaaatattcaatttattatttacaggcaattttaaactcaaaatatttagagtggtttagtgccttacatggagaagttttccgaggtggttatattgctaggggaactaaggtgcttaagaacttgcctattaggaaaattgattttgataatcttgaagaagcaaatctacatgatctaattgcgaccaagcaaaaagagcttatagagatttatgacaaaatagatgttaatgtaaataataaaagagttctgaccccattgcaacgtatgtttaaacgagagaaagaggttttagaccaattgttgagtcgactgtataacttaggtgtagatgattccttgatcccttatattaaggatttgtatgaagctcattaaGGATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAATAATTGCTCATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTGCCATCCAGCTGATATCCCCTATAGTGAGTCGTATTACATGGTCATAGCTGTTTCCTGGCAGCTCTGGCCCGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGCCTCCTCTAGACCAGCCAGGACAGAAATGCCTCGACTTCGCTGCTGCCCAAGGTTGCCGGGTGACGCACACCGTGGAAACGGATGAAGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGGGGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAAACATCATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATCGAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGCGGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAAACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAACATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAGGATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCTGGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGCAAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTATCAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCCTCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTAGTCGGCAAATAACCCTCGAGCCACCcatgaccaaaatcccttaacgGCATGCgcaccgccggacatcagcgctagcggagtgtatactggcttactatgttggcactgatgagggtgtcagtgaagtgcttcatgtggcaggagaaaaaaggctgcaccggtgcgtcagcagaatatgtgatacaggatatattccgcttcctcgctcactgactcgctacgctcggtcgttcgactgcggcgagcggaaatggcttacgaacggggcggagatttcctggaagatgccaggaagatacttaacagggaagtgagagggccgcggcaaagccgtttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaagaagatcatcttattaatcagataaaatatttctagAGGCCTcccctgattctgtggataaccGTcctaggTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctggaattcatgAtacgaa";

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

Application.Plasmid.directive('plasmidEditor', ['$parse', '$timeout', '$filter', '$compile', function($parse, $timeout, $filter, $compile) {

    return {
        restrict: "A",
        require:'ngModel',
        //rep1ace: true,
        link: function(scope, element, attrs, ngModel) {
            attrs.$set('spellcheck', false);

            scope.features = $parse(attrs.features)(scope) || [];

            //generate the mirror
            var parent = element.parent(),
                output = angular.element('<pre>'),
                label = angular.element('<label>');

            parent.append(output);
            parent.append(label);
            element[0].parentNode.replaceChild(parent, element);
            label.append(element);
            parent[0].className = 'ldt ' + element[0].className;

            /* key functions  */

            function genFiltered (text) {
                text = typeof text != 'undefined' ? text : ngModel.$modelValue;
                return $filter('features')(text, scope.features);
            }


            /* updating view */

            function setOutput (html) {
                output[0].innerHTML = html;
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
Application.Plasmid.filter('features', [function() {
    return function (text, features) {
        if (features && angular.isArray(features)) {
            var html = text;
            text = html.replace(/(<([^>]+)>)/ig, html);
            
            console.log(html);
            console.log(text);


            //future - pull tags out of HTML and save locations


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

                    angular.forEach(value, function(feat) {
                        var featName = angular.lowercase(feat.label).replace(/[ _]/gi, '');

                        //console.log(feat);
                        //console.log(key + " " + feat.label);

                        if (key == 'start') {
                            newText = newText.slice(0, index) +
                                '<span' + featName + ' ' +
                                'name="' + feat.label + '" ' +
                                'style="' +
                                (feat.css.color ? "color: " + feat.css.color + ";" : "") +
                                (feat.css.background ? "background-color: " + feat.css.background : "") +
                                '">' +
                                newText.slice(index);
                        } else {
                            newText = newText.slice(0, index) +
                                '</span' + featName + '>' +
                                newText.slice(index);
                        }
                    });

                })
            });

            console.log(newText);

            return newText;
        } else {
            return text;
        }
    };
}]);

Application.Plasmid.directive('feature', [function() {
    return {
        restrict : 'E',
        link: function(scope, element, attrs, controller) {
            //borrow from angularUI tooltip?
            element.className = element.className + ' feature';
        }
    }
}]);