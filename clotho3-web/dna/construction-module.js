"use strict";

Application.Dna.service('Construction', ['Clotho', 'DNA', 'Digest', 'PCR', '$parse', '$q', function(Clotho, DNA, Digest, PCR, $parse, $q) {

    //testing
    var dnaModules = {};
    dnaModules.DNA = DNA;
    dnaModules.Digest = Digest;
    dnaModules.PCR = PCR;

    var process = function(file) {
        //don't alter the construction file's dictionary
        var dict = file.dictionary;

        //process the dictionary, make all objects
        var dictPromises = {};
        angular.forEach(dict, function(value, key) {
            if (angular.isString(value)) {
                dict[key] = {computed : false, value : value};
            }

            if (!!value.Clotho) {
                //testing clientSide for now
                var parsed = $parse(value.value)(dnaModules);
                dictPromises[key] = {Clotho : true, computed : false, value : parsed};

                //future - dictPromises[key] = Clotho.submit(value.value);
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
                    
                    console.log('\n\n\n\ninput:', input);
                    
                    //future - handle 2+ layers deep
                    //for now, go one layer deep
                    if (!angular.isString(input)) {
                        //assume array or object
                        var parsedInput = (angular.isArray(input)) ? [] : {};
                        angular.forEach(input, function(child) {
                            //todo - handle object scenario
                            parsedInput.push(($parse(child)(dict)).value);
                        });
                        console.log(parsedInput);
                        inputs.push(parsedInput);
                    } else {
                        parsedInput = ($parse(input)(dict)).value;
                        console.log(parsedInput);
                        inputs.push(parsedInput);
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
                dict[step.output] = {computed : true, value : result};

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
            return dict.final.value
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
                '<tr ng-repeat="(key, value) in dict" ng-class="{\'info\' : value.computed}">' +
                '<td> <code ng-bind="key"></code> </td>' +
                '<td ng-bind="(value.value) | json"></td>' +
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
        template: '<div class="control-group span4" ng-class="{\'warning\' : hasInvalidItem }">' +
            '<label class="control-label">{{ fieldName }}</label>' +
            '<div class="controls">' +
            '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-change="checkInput()" disabled/>' +
            '</div>' +
            '</div>',
        link : function(scope, element, attrs) {
            scope.hasInvalidItem = false;
            scope.checkInput = function() {
                scope.hasInvalidItem = !!scope.dict[scope.model];
            };
        }
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

Application.Dna.directive('constructionStep', ['$parse', '$compile', '$http', function($parse, $compile, $http) {
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
                        '<div class="output">Output: <code>{{ step.output }}</code></div>' +
                        '</div>' +
                        '<div class="arrow">' +
                        '</div>');


                    function generateConstructionField (type, name, model) {
                        switch(angular.lowercase(type)) {
                            case 'value' : {
                                return '<div construction-field-string="'+name+'" ng-model="'+model+'"></div>';
                            }
                            case 'boolean' : {
                                return '<div construction-field-boolean="'+name+'" ng-model="'+model+'"></div>';
                            }
                            case 'array' : {
                                return '<div construction-field-array="'+name+'" ng-model="'+model+'" construction-dictionary="dict"></div>';
                            }
                            default : {
                                return '<div construction-field="'+name+'" ng-model="'+model+'" construction-dictionary="dict"></div>';
                            }
                        }
                    }

                    //if pass in scope.fields, assume want to overrun templates
                    if (angular.isDefined(scope.fields)) {
                        angular.forEach(scope.fields, function(field) {
                            field.model = 'step.' + field.model;

                            var html = generateConstructionField(field.type, field.name, field.model);
                            template.find('stepFields').append(html)
                        });
                    }
                    else {
                    //try to get the template for the reaction
                        $http.get('/dna/construction/steps/' + scope.step.reaction + '.html')
                        .then(function(data) {
                            //works - great
                            console.log(data.data);
                            console.log(scope);
                            var fieldsTemplate = $compile(data.data)(scope);
                            template.find('stepFields').append(fieldsTemplate);

                        }, function(data) {
                            //didn't  work - try to parse input fields
                            console.log('template not found');

                            angular.forEach(scope.step.input, function(input, index) {
                                var html;

                                if (angular.isArray(input)) {
                                    html = generateConstructionField('array', '', 'step.input['+index+']');
                                }
                                else if (input == "true" || input == "false") {
                                    html = generateConstructionField('boolean', '', 'step.input['+index+']');
                                }
                                else {
                                    html = generateConstructionField('value', '', 'step.input['+index+']');
                                }

                                template.find('stepFields').append(html)
                            });
                        });
                    }

                    element.html($compile(template)(scope))
                },
                post: function postLink(scope, element, attrs, ngModel) {

                    scope.showReaction = true;
                    scope.toggleReaction = function() {
                        scope.showReaction = !scope.showReaction;
                    }
                }
            }
        }
    }
}]);