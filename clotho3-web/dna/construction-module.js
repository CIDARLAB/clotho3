"use strict";

Application.Dna.service('Construction', ['Clotho', 'DNA', 'Digest', 'PCR', '$parse', '$q', function(Clotho, DNA, Digest, PCR, $parse, $q) {

    //testing
    var dnaModules = {};
    dnaModules.DNA = DNA;
    dnaModules.Digest = Digest;
    dnaModules.PCR = PCR;

    var process = function(file) {
        //don't alter the construction file's dictionary
        var dict = angular.copy(file.dictionary);

        //process the dictionary
        var dictPromises = {};
        angular.forEach(dict, function(value, key) {
            if (!!value.Clotho) {
                //dictPromises[key] = Clotho.submit(value.value);

                //testing clientSide for now
                dictPromises[key] = $parse(value.value)(dnaModules);
            }
        });


        var stepsChain = $q.defer();
        angular.forEach(file.steps, function(step) {

            //run function next in chain (so dict populated)
            stepsChain.promise.then(function() {

                //new deferred for each step -- or return something in Clotho.run block
                var stepDeferred = $q.defer();

                //parse inputs
                //future - handle promises in inputs
                //todo - handle Clotho specifics after parse (e.g. enzyme)
                var inputs = [];
                angular.forEach(step.input, function(input){
                    //future - handle 2+ layers deep
                    //for now, go one layer deep
                    if (!angular.isString(input)) {
                        //assume array or object
                        var parsedInput = (angular.isArray(input)) ? [] : {};
                        angular.forEach(input, function(child) {
                            //todo - handle object scenario
                            parsedInput.push($parse(child)(dict));
                        });
                        inputs.push(parsedInput);
                    } else {
                        inputs.push($parse(input)(dict));
                    }
                });

                console.log('running ' + step.reaction + ' with inputs:', inputs);

                //future -  run in clotho
                /*
                Clotho.run(step.reaction, inputs).then(function(result) {
                    dict[step.output] = result;
                    console.log(dict);
                    stepDeferred.resolve()
                });
                */

                //note - for now running client side
                //use apply to match format on server
                var result = $parse(step.reaction)(dnaModules).apply(null, inputs);
                console.log('result of '+step.reaction+':', result);
                dict[step.output] = result;

                return stepDeferred.promise;
            });
        });

        //kickoff
        $q.all(dictPromises).then(function(processedDict) {
            angular.extend(dict, processedDict);
            stepsChain.resolve();
        });

        //at the end
        return stepsChain.promise.then(function() {
            console.log('final: ' + dict.final, dict);
            return dict.final
        });


    };

    return {
        process : process
    }


}]);

Application.Dna.directive('constructionDictionaryView', [function() {
    return {
        restrict : 'EA',
        require : 'ngModel',
        scope : {
            dict : '=ngModel'
        },
        template : '<table class="table table-bordered table-striped table-condensed">' +
            '<thead><tr><th>Key</th> <th>Value</th></tr></thead>' +
            '<tbody>' +
                '<tr ng-repeat="(key, value) in dict">' +
                '<td> <code contenteditable ng-model="key"></code> </td>' +
                '<td contenteditable ng-model="value"></td>' +
                '</tr>' +
                '<tr><td colspan="2"><button class="btn" ng-click="addTerm()">Add new Key</button></td></tr>' +
            '</tbody>' +
            '</table>',
        link : function(scope, element, attrs, ngModel) {

            scope.$watch('dict', function(newval) {
                console.log(newval);
                ngModel.$render()
            }, true);

            scope.addTerm = function() {
                angular.extend(scope.dict, {"" : ""});
            }
        }
    }
}]);

Application.Dna.directive('constructionField', ['$parse', function($parse) {
    return {
        restrict : "EA",
        require: "ngModel",
        replace: true,
        scope: {
            fieldName : '@constructionField',
            model : '=ngModel',
            dict : '=constructionDictionary'
        },
        template: '<div class="control-group span4">' +
            '<label class="control-label">{{ fieldName }}</label>' +
                '<div class="controls">' +
                    '<select class="span12" placeholder="{{fieldName}}" ng-model="model" ng-options="key as key for (key,val) in dict" disabled></select>' +
                '</div>' +
            '</div>',
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {
                    /*
                    console.log(scope);
                    console.log('model before: ' + scope.model);
                    scope.model= $parse('step.' + scope.model)(scope);
                    console.log('model :' + scope.model);
                    */
                },
                post: function postLink(scope, element, attrs, ngModel) {
                }
            }
        }
    }
}]);

Application.Dna.directive('constructionFieldString', ['$parse', function($parse) {
    return {
        restrict : "EA",
        require: "ngModel",
        replace: true,
        scope: {
            fieldName : '@constructionFieldString',
            model : '=ngModel'
        },
        template: '<div class="control-group span4">' +
            '<label class="control-label">{{ fieldName }}</label>' +
            '<div class="controls">' +
            '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" disabled/>' +
            '</div>' +
            '</div>',
        link : function(scope, element, attrs) {}
    }
}]);

Application.Dna.directive('constructionFieldBoolean', ['$parse', function($parse) {
    return {
        restrict : "EA",
        require: "ngModel",
        replace: true,
        scope: {
            fieldName : '@constructionFieldBoolean',
            model : '=ngModel'
        },
        template: '<div class="control-group span4">' +
            '<label class="control-label">{{ fieldName }}</label>' +
            '<div class="controls">' +
            '<input class="span12" type="checkbox" ng-model="model" ng-true-value="true" ng-false-value="false" disabled/>' +
            '</div>' +
            '</div>',
        link : function(scope, element, attrs) {}
    }
}]);

Application.Dna.directive('constructionFieldArray', ['$parse', function($parse) {
    return {
        restrict : "EA",
        require: "ngModel",
        replace: true,
        scope: {
            fieldName : '@constructionFieldArray',
            model : '=ngModel',
            dict : '=constructionDictionary'
        },
        template: '<div class="control-group span4" ng-class="{\'error\' : hasInvalidItem }">' +
            '<label class="control-label">{{ fieldName }}</label>' +
            '<div class="controls">' +
            '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-list ng-change="checkInput()" disabled/>' +
            '</div>' +
            '</div>',
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {
                },
                post: function postLink(scope, element, attrs, ngModel) {
                    scope.hasInvalidItem = false;
                    scope.checkInput = function() {
                        var checkArrayRound = false;
                        angular.forEach(scope.model, function(item) {
                            if (!scope.dict[item]) {
                                checkArrayRound = true;
                            }
                        });
                        scope.hasInvalidItem = !!(checkArrayRound);
                    };
                }
            }
        }
    }
}]);

Application.Dna.directive('constructionStep', ['$parse', '$compile', function($parse, $compile) {
    return {
        restrict : "EA",
        require: "ngModel",
        scope: {
            step : '=ngModel',
            dict : '=constructionDictionary',
            fields : '=constructionStep'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {

                    var template = angular.element('<div class="constructionStep clearfix">'+
                        '<div class="reaction"><i ng-class="{\'icon-resize-full\' : !showReaction, \'icon-resize-small\' : showReaction }" ng-click="toggleReaction()"></i> {{ step.reaction }}</div>'+
                        '<div ng-show="showReaction">' +
                        //'{{fields}}' +
                        '<div class="steps">' +
                        '<stepFields class="row-fluid clearfix"></stepFields>' +
                        //'{{ step }}'+
                        '</div>' +
                        '<div class="output" tooltip="{{ dict[step.output] }}" tooltip-append-to-body="true">Output: <code>{{ step.output }}</code></div>' +
                        '</div>' +
                        '<div class="arrow">' +
                        '</div>');

                    angular.forEach(scope.fields, function(field) {
                        field.model = 'step.' + field.model;

                        console.log(field.type);

                        var html;

                        switch(angular.lowercase(field.type)) {
                            case 'value' : {
                                html = '<div construction-field-string="'+field.name+'" ng-model="'+field.model+'"></div>';
                            break;
                            }
                            case 'boolean' : {
                                html = '<div construction-field-boolean="'+field.name+'" ng-model="'+field.model+'"></div>';
                            break;
                            }
                            case 'array' : {
                                html = '<div construction-field-array="'+field.name+'" ng-model="'+field.model+'" construction-dictionary="dict"></div>';
                            break;
                            }
                            default : {
                                html = '<div construction-field="'+field.name+'" ng-model="'+field.model+'" construction-dictionary="dict"></div>';
                            }
                        }

                        template.find('stepFields').append(html)
                    });

                    element.html($compile(template)(scope))
                },
                post: function postLink(scope, element, attrs, ngModel) {

                    //extend dictionary with output term
                    var dictAddition = {};
                    dictAddition[scope.step.output] = {
                        "generated" : true,
                        "value" : scope.step.output
                    };
                    angular.extend(scope.dict, dictAddition);

                    scope.showReaction = true;

                    scope.toggleReaction = function() {
                        scope.showReaction = !scope.showReaction;
                    }
                }
            }
        }
    }
}]);