'use strict';

Application.Plasmid.controller('PlasmidCtrl', ['$scope', function($scope) {
    $scope.sequence = "ACAGTCACTGACTAGCTACGTACATCATCATCACAGCGTAGCTAGCTAGCTAGCATCGATCGATCGATCGAC";

    $scope.featureList = [
        {
            "pos" : {
                "start" : 6,
                end : 12
            },
            "css" : {
                "color" : "",
                "background" : "#99EEEE"
            }
        },
        {
            "match" : "ACAG",
            "css" : {
                "color" : "#BB99BB",
                "background" : ""
            }
        },
        {
            "match" : "TAGCTAG",
            "css" : {
                "color" : "",
                "background" : "#EECC99"
            }
        }

    ];

    $scope.emptyFeat = function() {
        return {pos: '', match: '', css : {color: '', background: '' }};
    };

    $scope.reg_match = /^[acgtACGT]+$/;
    $scope.reg_color = /^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/;
    $scope.reg_pos = /^[0-9]+\-[0-9]+$/;

    $scope.featureValid = function(feat) {
        return (
                (feat.match != '' && angular.isDefined(feat.match)) ?
                $scope.reg_match.test(feat.match) :
                (angular.isDefined(feat.pos) && $scope.reg_pos.test(feat.pos))
            ) &&
            (angular.isDefined(feat.css.color) ?
                ($scope.reg_color.test(feat.css.color)) :
                (angular.isDefined(feat.css.background) && ($scope.reg_color.test(feat.css.background))))
    };

    $scope.addFeature = function(feat) {
        if ($scope.featureValid(feat)) {
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
            //todo - force apply on plasmid-editor
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

            function genFiltered () {
                return $filter('features')(ngModel.$modelValue, scope.features);
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

//todo - don't manipulate directly - create map
Application.Plasmid.filter('features', [function() {
    return function (text, features) {
        if (features && angular.isArray(features)) {
            text = text.toString();
            angular.forEach(features, function(feat, ing) {
                //check if position
                if (feat.pos) {
                    var start = text.substring(0, feat.pos.start),
                        highlight = text.substring(feat.pos.start, feat.pos.end),
                        end = text.substring(feat.pos.end);

                    text = start +
                        '<span style="' +
                        (feat.css.color ? "color: " + feat.css.color + ";" : "") +
                        (feat.css.background ? "background-color: " + feat.css.background : "") +
                        '">'+highlight+'</span>' +
                        end;

                } else {
                    text = text.replace(new RegExp(feat.match, 'gi'), '<span style="' +
                        (feat.css.color ? "color: " + feat.css.color + ";" : "") +
                        (feat.css.background ? "background-color: " + feat.css.background : "") +
                        '">$&</span>');
                }
            });
            return text;
        } else {
            return text;
        }
    };
}]);