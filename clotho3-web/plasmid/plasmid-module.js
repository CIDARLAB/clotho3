'use strict';

//TODO - PLASMID PROVIDER -- validation, store features, add features, etc.

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

//future - make provider, not singleton service
Application.Plasmid.service('Plasmid', ['$window', function($window) {

    var features = [
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

    var regexps = {
        reg_match : /^[acgtACGT]+$/,
        reg_color : /^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/,
        reg_pos : /^[0-9]+\-[0-9]+$/
    };

    var emptyFeature = function() {
        return {label: "", pos: '', match: '', css : {color: '', background: '' }, active: true};
    };

    var randomColor = function() {
        return '#'+Math.floor(Math.random()*16777215).toString(16);
    };

    var featureValid = function(feat) {
        return (
            (feat.match != '' && angular.isDefined(feat.match)) ?
                regexps.reg_match.test(feat.match) :
                (angular.isDefined(feat.pos) && regexps.reg_pos.test(feat.pos))
            ) &&
            ( (angular.isDefined(feat.css.color) && feat.css.color != '') ?
                (regexps.reg_color.test(feat.css.color)) :
                (angular.isDefined(feat.css.background) && (regexps.reg_color.test(feat.css.background))))
    };

    var addFeature = function(feat) {
        if (featureValid(feat)) {

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
                feat.css.background = randomColor();

            console.log(feat);

            features.push(feat);
        } else {
            console.log("invalid");
        }
    };

    return {
        features : features,
        regexps : regexps,
        emptyFeature : emptyFeature,
        randomColor : randomColor,
        featureValid : featureValid,
        addFeature : addFeature
    }
}]);


Application.Plasmid.controller('PlasmidCtrl', ['$scope', '$window', '$document', 'Plasmid', function($scope, $window, $document, Plasmid) {
    $scope.sequence = "GATCTgttgacggctaGCTCAGTCCTAGGTagctACAGTGCTAGCTCTCTGGAGATTAACGAGGAGAAATACTAGATGGTTCATGATCATAAgcttgaattagccaaacttattcgcaactatgagacgaatagaaaagaatgtctaaattccagatataatgaaacacttttacgaagtgattatcttgatccattttttgaacttcttggctgggatattaaaaataaagctggaaaaccgactaatgaaagagaggttgtcttggaagaggcacttaaagcaagtgcatctgaacattctaaaaaaccagattatacattcagacttttttctgaaagaaagtttttcttggaagctaaaaaaccatcagttcatattgaatcggataatgaaactgctaaacaagtgcgaagatatggctttaccgccaaactaaaaatttcagttttatcaaattttgaatatttagttatttatgatacctctgtaaaggttgatggtgatgatacctttaataaggcacgtataaaaaaataccattacacagagtatgaaactcactttgatgaaatttgtgacttattaggaagagagtccgtttactctgggaattttgataaagaatggttgagtatcgaaaataaaattaatcacttttctgtagataccttatttttaaaacagattaatacatggcgtctattgcttggtgaagaaatctataagtatcaacctacgatacaagagaatgagcttaatgacattgtacagagctatctgaatagaattAGATCTatttttttgagagtctgtgaagatagaaatttagagacttatcagacattactgaattttgcttcaagtaatgatttctccgctcttattgataagtttaagcaggcagatcgttgctataattcaggcctatttgatcaattgcttacagagcaaattattgaggatattagttctgtattttgggtaatcattaagcaattatattatccagaaagtccttattcatttagtgtgttctcttcggatattttaggtaatatttacgaaatatttttatctgagaaattagtaattaatcaaagcagagttgagttagtcaagaaaccagagaatttagatagagacattgtcacaacaccaacctttattattaatgacatcttgagaaatacggttctaccgaagtgctatggaaaaacagatatagaaattctacagctaaaatttgctgatattgcttgtggttcgggagcatttttactggagttgttccaattacttaatgatactctagttgactattatttaagtagtgatacttctcaattaattccaacaggtatcggtacttataagctgtcttatgaaatcaagagaaaggttctattaagttgtatttttggcatagataaggacttaaatgctgtagaggctgcaaagttcggattgttgctaaaattattagagggtgaagacgtacaatctatagctaatattagaccagttctcccagatttattagataacatactttttggtaacagtttattagaaccagaaaaagtcgagcttgatcatcaggtagaagtaaatccgttagatttttcGGATCTGaaAtttgatgtaattgttggcaaccctccatatatgaaatcagaggatatgaagaatattactcctttggagttacctttatataagaaaaactatgtttctgcttataagcaatttgataaatatttcttgttcttagagcggggtttagctctattaaaagaagagggaatacttggatatattgttccaagtaaatttactaaagtgggtgcagggaaaaagttacgggaattactaacagataagggttatcttgactctattgtttcttttggtgctaatcaaatatttcaggataaaacaacttatacttgtttacttattttaagaaaaactccAcatactgattttaaatatgcagaggttcgtaatttaattgactggaaagtgcgtaaagctgatgctatggaattttcctctcaacaactgagtacattgcaaagtgatgcgtggattttaattccatctgaattaatctcagtttatcatcagatattagcacaaagccaaaagctagaggatattgtcggtattgataatatatttaatgggattcaaaccagtgctaatgatgtctatatttttgtgccaactcatgaggatactgaaaactattattttataaagaaaggacaagagtacaaaattgaaaaggaaattacgaagccttattttaaaacaacgagtggtgaggataacttatatacttaccgtactttcaagcctaatgcccgagtcatttatccgtatactcaaactgagagtagtgtagaactaattcctttagatgaaatacgagaaatttttcctttagcatacaaatatttaatgtcgcttaagttcgttttaagtagccccaaacgagatataaaacctagacctaaaacaacaaatgaatggcataggtatggacggcatcaaagtctcgataattgtgggttgagtcagaaaattattgtaggtgtgctttcagttggtgataagtacgctatagatacttatggaacgttgatttcatcaggcggtacggctggatactgtgtggttgctcttccagatgattgtaaatattcaatttattatttacaggcaattttaaactcaaaatatttagagtggtttagtgccttacatggagaagttttccgaggtggttatattgctaggggaactaaggtgcttaagaacttgcctattaggaaaattgattttgataatcttgaagaagcaaatctacatgatctaattgcgaccaagcaaaaagagcttatagagatttatgacaaaatagatgttaatgtaaataataaaagagttctgaccccattgcaacgtatgtttaaacgagagaaagaggttttagaccaattgttgagtcgactgtataacttaggtgtagatgattccttgatcccttatattaaggatttgtatgaagctcattaaGGATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAATAATTGCTCATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTGCCATCCAGCTGATATCCCCTATAGTGAGTCGTATTACATGGTCATAGCTGTTTCCTGGCAGCTCTGGCCCGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGCCTCCTCTAGACCAGCCAGGACAGAAATGCCTCGACTTCGCTGCTGCCCAAGGTTGCCGGGTGACGCACACCGTGGAAACGGATGAAGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGGGGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAAACATCATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATCGAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGCGGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAAACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAACATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAGGATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCTGGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGCAAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTATCAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCCTCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTAGTCGGCAAATAACCCTCGAGCCACCcatgaccaaaatcccttaacgGCATGCgcaccgccggacatcagcgctagcggagtgtatactggcttactatgttggcactgatgagggtgtcagtgaagtgcttcatgtggcaggagaaaaaaggctgcaccggtgcgtcagcagaatatgtgatacaggatatattccgcttcctcgctcactgactcgctacgctcggtcgttcgactgcggcgagcggaaatggcttacgaacggggcggagatttcctggaagatgccaggaagatacttaacagggaagtgagagggccgcggcaaagccgtttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaagaagatcatcttattaatcagataaaatatttctagAGGCCTcccctgattctgtggataaccGTcctaggTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctggaattcatgAtacgaa";

    $scope.editable = true;
    $scope.textSelected = false;

    $scope.features = Plasmid.features;

    angular.forEach($scope.features, function(feat) {
        feat.active = !!feat.active || true;
    });

    //todo - expose, make editable
    $scope.reg_match = Plasmid.regexps.reg_match;
    $scope.reg_color = Plasmid.regexps.reg_color;
    $scope.reg_pos = Plasmid.regexps.reg_pos;


    $scope.emptyFeat = Plasmid.emptyFeature;

    $scope.logFeatures = function() {console.log($scope.features)};

    $scope.featureValid = Plasmid.featureValid;

    $scope.toggleActive = function(feat, val) {

        if (angular.isArray(feat)) {
            angular.forEach(feat, function(cur) {
                cur.active = val || !cur.active;
            });
        } else {
            feat.active = val || !feat.active;
        }
        $scope.$broadcast('Plasmid:Process');
    };

    $scope.randomColor = Plasmid.randomColor;

    $scope.addFeature = function(feat) {
        Plasmid.addFeature(feat);
        $scope.new = Plasmid.emptyFeature();
    };

}]);

Application.Plasmid.directive('plasmidEditor', ['$parse', '$timeout', '$filter', '$compile', '$document', '$window', 'Plasmid', function($parse, $timeout, $filter, $compile, $document, $window, Plasmid) {

    return {
        restrict: "A",
        require:'ngModel',
        scope : {
            features: '=',
            editable: '@',
            sequence: '=ngModel'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, controller) {
                    attrs.$set('spellcheck', false);
                    attrs.$set('contenteditable', true); //todo - check controller

                    element.css({
                        'font-size': '16px',
                        'font-family': 'monospace',
                        'color': 'black',
                        'outline': 'none',
                        'resize': 'none',
                        'min-width': '100%',
                        'overflow': 'hidden',
                        'box-sizing': 'border-box',
                        '-moz-padding-start': '1px',
                        'min-height': '22px',
                        'position': 'relative'
                    });

                    element.parent().prepend($compile(angular.element('<div class="btn-group"><button class="btn btn-small" ng-click="highlight()" ng-class="{\'btn-info\' : ngModel.$dirty}" ng-disabled="ngModel.$pristine">Process</button><button class="btn btn-small" ng-click="addFeatureSelection()" ng-disabled="!textSelected">Annotate Selection</button></div>'))(scope));
                },
                post: function(scope, element, attrs, ngModel) {
                    /* key functions  */

                    //for $pristine etc. access in template
                    scope.ngModel = ngModel;

                    //console.log(scope);

                    function genFiltered (text) {
                        text = typeof text != 'undefined' ? text : ngModel.$viewValue;
                        return $filter('features')(text, scope.features);
                    }

                    /* updating view */

                    function setOutput (html) {
                        element.html(html);
                        $compile(element.contents())(scope);
                        ngModel.$setPristine();
                    }

                    scope.highlight = function() {
                        setOutput(genFiltered());
                    };

                    scope.logModel = function() {
                        console.log(ngModel.$viewValue);
                        console.log(ngModel.$modelValue);
                    };

                    //todo - checks:
                    // select opp direction
                    // pos is for whole string (across text nodes)
                    scope.addFeatureSelection = function () {
                        var feature = Plasmid.emptyFeature();

                        if ($window.getSelection) {
                            var sel = $window.getSelection();
                            if (typeof sel != 'undefined' && sel.rangeCount) {

                                var container = angular.element("<div>");
                                for (var i = 0, len = sel.rangeCount; i < len; ++i) {
                                    container.append(sel.getRangeAt(i).cloneContents());
                                }
                                console.log(container);

                                //todo - see http://stackoverflow.com/questions/3597116/insert-html-after-a-selection

/*

                                var newEl = document.createElement("feat");
                                newEl.setAttribute('index', scope.features.length.toString());
                                sel.getRangeAt(0).surroundContents(newEl);
*/
                                //feature = container[0].innerHTML;

                                console.log(sel);






                                feature.pos = {"start" : sel.baseOffset, "end" : sel.extentOffset};

                                feature.css.background = '#'+Math.floor(Math.random()*16777215).toString(16);

                            }
                        }
                        else if (typeof $document.selection != "undefined") { //IE -- todo doesn't really work
                            if ($document.selection.type == "Text") {
                                feature = $document.selection.createRange().htmlText;
                            }
                        }
                        console.log(feature);
                        scope.features.push(feature);
                    };


                    //model -> view

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

                    scope.$on('Plasmid:Process', function() {
                        setOutput(genFiltered());
                    });


                    /* watchers */
                    //shallow watch - only with count
                    scope.$watchCollection('features', function(newVal, oldVal) {
                        if (!!newVal && !!oldVal) {
                            setOutput(genFiltered());
                        }
                    });

                    element.bind( "keyup change" , function() {
                        console.log('change to view model');
                        scope.$apply(
                            ngModel.$setViewValue(element.html())
                        );
                    });

                    scope.textSelected = false;
                    $document.bind('mouseup', function() {
                        if (typeof $window.getSelection != "undefined") {
                            scope.$apply(scope.textSelected = !!$window.getSelection().toString());
                        } else if (typeof $document.selection != "undefined" && $document.selection.type == "Text") {
                            scope.$apply(scope.textSelected = !!$document.selection.createRange().text);
                        }
                    });
                }
            }
        }
    }
}]);

Application.Plasmid.filter('features', [function() {
    return function (text, features) {
        if (features && angular.isArray(features) && angular.isString(text)) {

            //note - at this point, we have well-formed HTML (or text on first pass)
            var raw = text;
            text = raw.replace(/(<([^>]+)>)/ig, "");

            //console.log(raw);


            /*
            //todo - pull tags out of HTML and save locations
            //future - don't use jquery (for speed)
            var startStops = $(raw);
            startStops.find(':not(end, start)').children().unwrap();
            console.log(startStops);

            var test = startStops.not('start, end').each(function() {
                console.log(this);
                var content = $(this).text();
                $(this).replaceWith(content);
            });
            //startStops.find('span').contents().unwrap();

            console.log(startStops);
            console.log(test.html());

            console.log(angular.element(raw).find('start'));
            */

            //todo - rewrite so can use form <feat index="">

            //note - still contain ng-scope in attrs
            var startEnds = raw.replace(/\s*<(\/?)(\w+)([^>]*?)>\s*/ig,  function(j,b,a,c){
                return   ({start:1, end:1}[a]) ?   ("<"+b+a+c+">")  : "";
            });

            //console.log(startEnds);


            var overlap;
            var findStartEnd = /<(\w+) feat="(\d+)"[^>]*><\/\1>/ig;
            while ((overlap = findStartEnd.exec(startEnds)) != null) {
                //console.log(overlap);

                var ind = overlap.index,
                    type = overlap[1],
                    featIndex = overlap[2];

                var feature = features[featIndex];
                if (feature.pos) {
                    feature.pos[type] = ind;
                }

                startEnds = startEnds.replace(overlap[0], '');
                findStartEnd.lastIndex = 0;
            }
            //console.log(startEnds);






            //create location map
            var locations = {};
            angular.forEach(features, function(feat, featIndex) {
                if (feat.active === false) { //skip
                } else {
                    if (feat.match) {

                        //todo - check reverse direction too - note rev. direction

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
                }
            });
            //console.log(locations);


            //loop from end, saving tag locations
            var reversed = [];
            angular.forEach(locations, function(key, val) {
                reversed.unshift(val);
            });

            var newText = text,
                backlog = [];

            angular.forEach(reversed, function(index) {
                angular.forEach(locations[index], function(value, key) {
                    if (text.length < index) {
                        console.log("string too short");
                        return;
                    }

                    if (!value.length) {
                        return;
                    }
                    //console.log(index, locations[index], value);

                    angular.forEach(value, function(featIndex) {

                        //console.log(index, key, value, backlog);

                        if (key == 'start') {

                            var indices = (backlog.length > 1) ? backlog.join("-") : featIndex;
                            //console.log(indices);

                            newText = newText.slice(0, index) +
                                ((backlog.length > 1) ? '</annotation>' : '') +
                                '<start feat="' + featIndex + '"></start>' +
                                '<annotation index="'+ indices + '">' +
                                newText.slice(index);


                            //splice out featIndex
                            var ind = backlog.indexOf(featIndex);
                            backlog.splice(ind, 1);

                        } else {
                            backlog.push(featIndex);

                            newText = newText.slice(0, index) +
                                '</annotation>' +
                                '<end feat="' + featIndex + '"></end>' +
                                ((backlog.length > 1) ? '<annotation index="' + backlog.slice(0,-1).join("-") + '">' : '') +
                                newText.slice(index);
                        }
                        //console.log(newText);
                    });

                })
            });
            //console.log(newText);

            return newText;
        } else {
            return text;
        }
    };
}]);

Application.Plasmid.directive('annotation', ['$tooltip', function($tooltip) {

    return {
        restrict : 'EA',
        replace: false,
        scope: {
            index: '@'
        },
        transclude:true,
        template: '<span tooltip="{{ feature.label }}" tooltip-placement="mouse" tooltip-append-to-body="true" ng-transclude></span>',
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, controller) {
                    scope.features = scope.$parent.features;
                    var matches = scope.index.split('-');
                    scope.feature = angular.copy(scope.features[matches.pop()]);

                    if (matches.length) {
                        angular.extend(scope.feature.css, {'border-bottom' : '2px dotted red'});
                        for (var ind = 0; ind < matches.length; ind++) {
                            angular.extend(scope.feature.css, scope.features[matches[ind]].css);
                            scope.feature.label += ", " + scope.features[matches[ind]].label;
                        }
                    }

                },
                post: function(scope, element, attrs, ctrl) {
                    //attrs.$set('style', "background-color: #FF0000");

                    element.css(scope.feature.css);

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
        }
    }
}]);