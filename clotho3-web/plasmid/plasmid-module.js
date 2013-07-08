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
Application.Plasmid.service('Plasmid', ['$window', '$timeout', function($window, $timeout) {

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
            "label" : "Restriction Site",
            "match" : "TAGCTA",
            "css" : {
                "color" : "#FF88CC",
                "background" : ""
            }
        },
        {
            "label" : "My Feature",
            "match" : "TAGCTAG",
            "css" : {
                "color" : "",
                "background" : "#EECC99"
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
        //truly random, but probably want light colors only
        //return '#'+Math.floor(Math.random()*16777215).toString(16);

        var r = (Math.round(Math.random()* 127) + 127).toString(16);
        var g = (Math.round(Math.random()* 127) + 127).toString(16);
        var b = (Math.round(Math.random()* 127) + 127).toString(16);
        return '#' + r + g + b;
    };

    var featureValid = function(feat) {

        if (angular.isUndefined(feat)) return false;

        if (feat.match && feat.match != '' && angular.isDefined(feat.match)) {
            return regexps.reg_match.test(feat.match)
        } else {
            return angular.isDefined(feat.pos) && (regexps.reg_pos.test(feat.pos) || (angular.isNumber(feat.pos.start) && angular.isNumber(feat.pos.end)))
        }

        return false;
    };

    var addFeature = function(feat, skipCheck) {
        if (!skipCheck) {
            if (!featureValid(feat)) {
                console.log('invalid');
                return false;
            }
        }
        //console.log('adding');

        if (!feat.label) {
            var milli = new Date;
            milli = milli.getMilliseconds();
            feat.label = "Feature " + milli;
        }

        if (feat.match) {
            feat.match = angular.uppercase(feat.match);
            delete feat['pos'];
        }
        else {
            if (/^([0-9]+)\-([0-9]+)$/.test(feat.pos)) {
                var positions = feat.pos.match(/^([0-9]+)\-([0-9]+)$/);
                feat.pos = {};
                feat.pos.start = positions[1];
                feat.pos.end = positions[2];
            }

        }

        //if color defined, make sure valid
        if (feat.css.color != '') {
            if (!regexps.reg_color.test(feat.css.color)) feat.css.color = randomColor();
        } else if (feat.css.background != '') {
            if (!regexps.reg_color.test(feat.css.background)) feat.css.background = randomColor();
        }

        if (feat.css.color == '' && feat.css.background == '')
            feat.css.background = randomColor();

        features.push(feat);
    };

    var addSelection =  function (ngModel) {

        var feature = emptyFeature();

        var range, start, end, expandedSelRange;

        if ($window.getSelection) {
            var sel = $window.getSelection();
            if (typeof sel != 'undefined' && sel.rangeCount) {

                range = sel.getRangeAt(0);
                start = range.cloneRange();
                end = range.cloneRange();
                start.collapse(true);
                end.collapse(false);
                expandedSelRange = range.cloneRange();


                // Range.createContextualFragment() would be useful here but is
                // non-standard and not supported in all browsers (IE9, for one)
                // start
                var elstart = document.createElement("div");
                elstart.innerHTML = '<start feat="'+ features.length +'"></start>';
                var fragstart = document.createDocumentFragment(), node, lastNode;
                while ( (node = elstart.firstChild) ) {
                    lastNode = fragstart.appendChild(node);
                }
                start.insertNode(fragstart);

                //end
                var elend = document.createElement("div");
                elend.innerHTML = '<end feat="'+ features.length +'"></end>';
                var fragend = document.createDocumentFragment(), node, lastNode;
                while ( (node = elend.firstChild) ) {
                    lastNode = fragend.appendChild(node);
                }
                end.insertNode(fragend);

                // Preserve the selection
                if (lastNode) {
                    expandedSelRange.setEndAfter(lastNode);
                    sel.removeAllRanges();
                    sel.addRange(expandedSelRange);
                }
            }
        }
        //IE -- verify working?
        else if (typeof document.selection != "undefined") {

            range = document.selection.createRange();
            start = range.duplicate();
            end = range.duplicate();
            expandedSelRange = range.duplicate();

            start.collapse(true);
            end.collapse(false);

            start.pasteHTML('<start feat="'+features.length+'"></start>');
            end.pasteHTML('<end feat="'+features.length+'"></end>');

            expandedSelRange.setEndPoint("EndToEnd", range);
            expandedSelRange.select();
        }

        console.log(feature);

        //todo - this is hacky
        ngModel.$setViewValue($('[plasmid-editor]').html());
        feature.pos = {'start' : '', 'end' : ''};

        addFeature(feature, true);
    };

    return {
        features : features,
        regexps : regexps,
        emptyFeature : emptyFeature,
        randomColor : randomColor,
        featureValid : featureValid,
        addFeature : addFeature,
        addSelection : addSelection
    }
}]);


Application.Plasmid.controller('PlasmidCtrl', ['$scope', 'Clotho', '$window', '$document', 'Plasmid', '$route', function($scope, Clotho, $window, $document, Plasmid, $route) {

    //feature set up
    function activateFeatures() {
        angular.forEach($scope.features, function(feat) {
            feat.active = !!feat.active || true;
        });
    }

    if ($route.current.params.id) {
        $scope.id = $route.current.params.id;
        Clotho.get($scope.id).then(function(data) {
            $scope.NucSeq = data;
            //demo
            $scope.features = Plasmid.features;
            activateFeatures();
            $scope.$broadcast('Plasmid:Process');
        });
    } else {
        $scope.NucSeq = {
            "isRNA": false,
            "lastModified": {
                "$date": "2013-07-03T20:15:50.122Z"
            },
            "isDegenerate": false,
            "lowerArray": [
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false
            ],
            "isCircular": false,
            "schema_id": "org.clothocad.model.NucSeq",
            "isDeleted": false,
            "isLocked": false,
            "id": "51d48676300484c59c1bbe3f",
            "isSingleStranded": false,
            "sequence": "GATCTgttgacggctaGCTCAGTCCTAGGTagctACAGTGCTAGCTCTCTGGAGATTAACGAGGAGAAATACTAGATGGTTCATGATCATAAgcttgaattagccaaacttattcgcaactatgagacgaatagaaaagaatgtctaaattccagatataatgaaacacttttacgaagtgattatcttgatccattttttgaacttcttggctgggatattaaaaataaagctggaaaaccgactaatgaaagagaggttgtcttggaagaggcacttaaagcaagtgcatctgaacattctaaaaaaccagattatacattcagacttttttctgaaagaaagtttttcttggaagctaaaaaaccatcagttcatattgaatcggataatgaaactgctaaacaagtgcgaagatatggctttaccgccaaactaaaaatttcagttttatcaaattttgaatatttagttatttatgatacctctgtaaaggttgatggtgatgatacctttaataaggcacgtataaaaaaataccattacacagagtatgaaactcactttgatgaaatttgtgacttattaggaagagagtccgtttactctgggaattttgataaagaatggttgagtatcgaaaataaaattaatcacttttctgtagataccttatttttaaaacagattaatacatggcgtctattgcttggtgaagaaatctataagtatcaacctacgatacaagagaatgagcttaatgacattgtacagagctatctgaatagaattAGATCTatttttttgagagtctgtgaagatagaaatttagagacttatcagacattactgaattttgcttcaagtaatgatttctccgctcttattgataagtttaagcaggcagatcgttgctataattcaggcctatttgatcaattgcttacagagcaaattattgaggatattagttctgtattttgggtaatcattaagcaattatattatccagaaagtccttattcatttagtgtgttctcttcggatattttaggtaatatttacgaaatatttttatctgagaaattagtaattaatcaaagcagagttgagttagtcaagaaaccagagaatttagatagagacattgtcacaacaccaacctttattattaatgacatcttgagaaatacggttctaccgaagtgctatggaaaaacagatatagaaattctacagctaaaatttgctgatattgcttgtggttcgggagcatttttactggagttgttccaattacttaatgatactctagttgactattatttaagtagtgatacttctcaattaattccaacaggtatcggtacttataagctgtcttatgaaatcaagagaaaggttctattaagttgtatttttggcatagataaggacttaaatgctgtagaggctgcaaagttcggattgttgctaaaattattagagggtgaagacgtacaatctatagctaatattagaccagttctcccagatttattagataacatactttttggtaacagtttattagaaccagaaaaagtcgagcttgatcatcaggtagaagtaaatccgttagatttttcGGATCTGaaAtttgatgtaattgttggcaaccctccatatatgaaatcagaggatatgaagaatattactcctttggagttacctttatataagaaaaactatgtttctgcttataagcaatttgataaatatttcttgttcttagagcggggtttagctctattaaaagaagagggaatacttggatatattgttccaagtaaatttactaaagtgggtgcagggaaaaagttacgggaattactaacagataagggttatcttgactctattgtttcttttggtgctaatcaaatatttcaggataaaacaacttatacttgtttacttattttaagaaaaactccAcatactgattttaaatatgcagaggttcgtaatttaattgactggaaagtgcgtaaagctgatgctatggaattttcctctcaacaactgagtacattgcaaagtgatgcgtggattttaattccatctgaattaatctcagtttatcatcagatattagcacaaagccaaaagctagaggatattgtcggtattgataatatatttaatgggattcaaaccagtgctaatgatgtctatatttttgtgccaactcatgaggatactgaaaactattattttataaagaaaggacaagagtacaaaattgaaaaggaaattacgaagccttattttaaaacaacgagtggtgaggataacttatatacttaccgtactttcaagcctaatgcccgagtcatttatccgtatactcaaactgagagtagtgtagaactaattcctttagatgaaatacgagaaatttttcctttagcatacaaatatttaatgtcgcttaagttcgttttaagtagccccaaacgagatataaaacctagacctaaaacaacaaatgaatggcataggtatggacggcatcaaagtctcgataattgtgggttgagtcagaaaattattgtaggtgtgctttcagttggtgataagtacgctatagatacttatggaacgttgatttcatcaggcggtacggctggatactgtgtggttgctcttccagatgattgtaaatattcaatttattatttacaggcaattttaaactcaaaatatttagagtggtttagtgccttacatggagaagttttccgaggtggttatattgctaggggaactaaggtgcttaagaacttgcctattaggaaaattgattttgataatcttgaagaagcaaatctacatgatctaattgcgaccaagcaaaaagagcttatagagatttatgacaaaatagatgttaatgtaaataataaaagagttctgaccccattgcaacgtatgtttaaacgagagaaagaggttttagaccaattgttgagtcgactgtataacttaggtgtagatgattccttgatcccttatattaaggatttgtatgaagctcattaaGGATCCtaaCTCGAcgtgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatCAATTCGACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAATAATTGCTCATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTGCCATCCAGCTGATATCCCCTATAGTGAGTCGTATTACATGGTCATAGCTGTTTCCTGGCAGCTCTGGCCCGTGTCTCAAAATCTCTGATGTTACATTGCACAAGATAAAAATATATCATCATGCCTCCTCTAGACCAGCCAGGACAGAAATGCCTCGACTTCGCTGCTGCCCAAGGTTGCCGGGTGACGCACACCGTGGAAACGGATGAAGGCACGAACCCAGTGGACATAAGCCTGTTCGGTTCGTAAGCTGTAATGCAAGTAGCGTATGCGCTCACGCAACTGGTCCAGAACCTTGACCGAACGCAGCGGTGGTAACGGCGCAGTGGCGGTTTTCATGGCTTGTTATGACTGTTTTTTTGGGGTACAGTCTATGCCTCGGGCATCCAAGCAGCAAGCGCGTTACGCCGTGGGTCGATGTTTGATGTTATGGAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAAACATCATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATCGAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGCGGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAAACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAACATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAGGATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCTGGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGCAAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTATCAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCCTCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTAGTCGGCAAATAACCCTCGAGCCACCcatgaccaaaatcccttaacgGCATGCgcaccgccggacatcagcgctagcggagtgtatactggcttactatgttggcactgatgagggtgtcagtgaagtgcttcatgtggcaggagaaaaaaggctgcaccggtgcgtcagcagaatatgtgatacaggatatattccgcttcctcgctcactgactcgctacgctcggtcgttcgactgcggcgagcggaaatggcttacgaacggggcggagatttcctggaagatgccaggaagatacttaacagggaagtgagagggccgcggcaaagccgtttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaagaagatcatcttattaatcagataaaatatttctagAGGCCTcccctgattctgtggataaccGTcctaggTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCTCCGAATTGgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctggaattcatgAtacgaa",
            "description": "temporary description",
            "name": "nucseq",
            "isLinear": false,
            "className": "org.clothocad.model.NucSeq",
            "uuid": "51d48676300484c59c1bbe3f"
        };

        $scope.features = Plasmid.features;
        activateFeatures();
    }


    /* functionality */

    //todo - expose, make editable
    $scope.reg_match = Plasmid.regexps.reg_match;
    $scope.reg_color = Plasmid.regexps.reg_color;
    $scope.reg_pos = Plasmid.regexps.reg_pos;

    $scope.emptyFeat = Plasmid.emptyFeature;
    $scope.featureValid = Plasmid.featureValid;
    $scope.randomColor = Plasmid.randomColor;

    $scope.logFeatures = function() {console.log($scope.features)};

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

    $scope.cancelAdd = function() {
        $scope.addForm = false;
        $scope.editForm = false;
        $scope.new = Plasmid.emptyFeature();
    };

    //todo - move to service
    $scope.positionObjectToString = function(feat) {
        if (typeof feat.pos != 'undefined')
            feat.pos = feat.pos.start + "-" + feat.pos.end;
    };
    //todo - move to service
    $scope.positionStringToObject = function(feat) {
        if (typeof feat.pos != 'undefined') {
            var positions = feat.pos.match(/^([0-9]+)\-([0-9]+)$/);
            feat.pos = {};
            feat.pos.start = positions[1];
            feat.pos.end = positions[2];
        }
    };

    $scope.cancelEdit = function(feat) {
        $scope.positionStringToObject(feat);

        $scope.editForm = false;
        $scope.addForm = false;
        $scope.editFeat = undefined;
        $scope.editFeatInd = undefined;
    };

    $scope.showAddForm =function() {
        $scope.editForm = false;
        $scope.addForm = true;
    };

    $scope.editFeature = function(feat, index) {
        $scope.editForm = true;
        $scope.addForm = false;
        $scope.editFeat = feat;

        //todo - hacky - avoid index
        $scope.editFeatInd = index;
        $scope.positionObjectToString($scope.editFeat);
    };

    $scope.saveEditFeature = function(feat) {

        $scope.cancelEdit(feat);
        feat.edited = true;

        Plasmid.features[$scope.editFeatInd] = feat;

        $scope.$broadcast('Plasmid:Process');
    };

    $scope.addFeature = function(feat) {
        Plasmid.addFeature(feat);
        $scope.cancelAdd();
    };

}]);

Application.Plasmid.directive('plasmidEditor', ['$parse', '$timeout', '$filter', '$compile', '$document', '$window', 'Plasmid', function($parse, $timeout, $filter, $compile, $document, $window, Plasmid) {

    //todo - seems that need to pass in ngModel="sequence" if don't want to manually set within link.

    return {
        restrict: "A",
        require:'ngModel',
        scope : {
            features: '=',
            editable: '=',
            sequence: '=ngModel'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, controller) {
                    attrs.$set('spellcheck', false);
                    attrs.$set('contenteditable', scope.editable);

                    element.css({
                        'font-size': '16px',
                        'font-family': 'monospace',
                        'color': 'black',
                        'outline': 'none',
                        'resize': 'none',
                        'min-width': '100%',
                        'overflow': 'hidden',
                        'box-sizing': 'border-box',
                        'min-height': '22px',
                        'position': 'relative'
                    });

                    element.parent().prepend($compile(angular.element('<div class="btn-group"><button class="btn btn-small" ng-click="highlight()" ng-class="{\'btn-info\' : ngModel.$dirty}" ng-disabled="ngModel.$pristine">Process</button><button class="btn btn-small" ng-click="addSelection()" ng-disabled="!textSelected">Annotate Selection</button>'))(scope));

                    scope.textSelected = false;

                },
                post: function(scope, element, attrs, ngModel) {
                    /* key functions  */

                    //for $pristine etc. access in template
                    scope.ngModel = ngModel;

                    function genFiltered (text) {
                        //console.log(ngModel.$viewValue);
                        text = typeof text != 'undefined' ? text : ngModel.$viewValue;
                        return $filter('features')(text, scope.features);
                    }

                    /* updating view */

                    function setOutput (html) {
                        element.html(html);
                        $compile(element.contents())(scope);
                        ngModel.$setViewValue(element.html());
                        ngModel.$setPristine();
                    }

                    scope.highlight = function() {
                        setOutput(genFiltered());
                    };

                    scope.logModel = function() {
                        console.log(ngModel.$viewValue);
                        console.log(ngModel.$modelValue);
                    };

                    scope.addSelection = function() {
                        //scope.logModel();
                        Plasmid.addSelection(ngModel);
                        //scope.logModel();
                    };


                    /* updating / rendering */

                    //verify - only called at instantiation?
                    ngModel.$render = function() {
                        //console.log(ngModel.$viewValue);
                        console.log('render', scope.sequence, ngModel);
                        element.html(genFiltered(ngModel.$viewValue));
                    };

                    //todo - hacky
                    scope.$watch(function() {
                        return ngModel.$viewValue;
                    }, function(viewVal) {
                        if (!!scope.sequence && angular.isUndefined(ngModel.$viewValue) || (!!scope.sequence && ngModel.$viewValue == ''))
                            ngModel.$setViewValue(scope.sequence);

                        console.log('view watch', scope.sequence, ngModel);

                        //element.html(genFiltered(viewVal));
                    });

                    scope.$watch(function() {
                        return ngModel.$modelValue;
                    }, function(modelVal) {

                        console.log('model watch', scope.sequence, ngModel);

                        //element.html(genFiltered(viewVal));
                    });


                    /* watchers */
                    scope.$on('Plasmid:Process', function() {
                        setOutput(genFiltered());
                    });

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

                    $document.bind('mouseup', function() {
                        //todo - limit to element
                        if (typeof $window.getSelection != "undefined") {
                            scope.$apply(scope.textSelected = !!$window.getSelection().toString() && element.has($window.getSelection().anchorNode).length);
                        } else if (typeof $document.selection != "undefined" && $document.selection.type == "Text") {
                            scope.$apply(scope.textSelected = !!$document.selection.createRange().text && element.has($window.getSelection().anchorNode).length);
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

            //note - still contain ng-scope in attrs
            var startEnds = raw.replace(/\s*<(\/?)(\w+)([^>]*?)>\s*/ig,  function(j,b,a,c){
                return   ({start:1, end:1}[a]) ?   ("<"+b+a+c+">")  : "";
            });

            //console.log(startEnds);


            var overlap;
            var findStartEnd = /<(\w+) feat="(\d+)"[^>]*?><\/\1>/ig;
            while ((overlap = findStartEnd.exec(startEnds)) != null) {
                //console.log(overlap);

                var ind = overlap.index,
                    type = overlap[1],
                    featIndex = overlap[2];

                var feature = features[featIndex];
                if (feature.pos) {
                    if (feature.edited) {
                        feature.edited = false;
                    } else {
                        feature.pos[type] = ind;
                    }
                }

                startEnds = startEnds.replace(overlap[0], '');
                findStartEnd.lastIndex = 0;
            }
            //console.log(startEnds);

            //console.log(features);





            //create location map
            var locations = {};
            angular.forEach(features, function(feat, featIndex) {
                if (feat.hasOwnProperty("active") && feat.active === false) { //skip
                } else {
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

                        //check reverse direction too....

                        var reverseMatch = angular.lowercase(feat.match);
                        reverseMatch = reverseMatch.replace(/[actg]/g, function (m) {
                            return {
                                'a': 't',
                                'c': 'g',
                                'g': 'c',
                                't': 'a'
                            }[m];
                        });
                        //if reverse complementary, don't check again

                        console.log(reverseMatch, feat.match);
                        
                        if (reverseMatch == angular.lowercase(feat.match).split("").reverse().join("")){
                        } else {
                            var reverseSeq = (angular.lowercase(text)).split("").reverse().join("");

                            for (var index, offset = 0, search = reverseSeq;
                                 (index = search.indexOf(reverseMatch, offset)) > -1;
                                 offset = index + feat.match.length
                                ) {
                                //todo - better logic?
                                //start
                                var end = reverseSeq.length - index;
                                locations[end] ?
                                    locations[end]['end'].push(featIndex) :
                                    locations[end] = {"start" : [], "end" : [featIndex]};
                                //end
                                var start = end - feat.match.length;
                                locations[start] ?
                                    locations[start]['start'].push(featIndex) :
                                    locations[start] = {"start" : [featIndex], "end" : []};

                            }
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