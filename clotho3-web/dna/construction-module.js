"use strict";

Application.Dna.service('Construction', ['Clotho', 'DNA', 'Digest', 'PCR', '$parse', '$q', function(Clotho, DNA, Digest, PCR, $parse, $q) {

    //testing
    var dnaModules = {};
    dnaModules.DNA = DNA;
    dnaModules.Digest = Digest;
    dnaModules.PCR = PCR;

    var process = function(file, returnDictionary) {
        //don't alter the construction file's dictionary
        var dict = file.dictionary;

        console.log('PROCESSING CONSTRUCTION FILE', file);

        dict = _.pick(dict, function(value, key) {
            return !value.computed
        });

        //console.log(file);
        //console.log(dict);
        //console.log(_.keys(dict));


        //process the dictionary, make all objects
        var dictPromises = {};
        angular.forEach(dict, function(value, key) {
            if (angular.isString(value)) {
                dict[key] = {key : key, computed : false, value : value};
            }

            if (!!value.Clotho && !(value.retrieved === true)) {
                //testing clientSide for now
                var parsed = $parse(value.value)(dnaModules);
                dictPromises[key] = {key: key, Clotho : true, retrieved : true, computed : false, value : parsed};

                //future - dictPromises[key] = Clotho.submit(value.value);
            }
        });


        //kickoff
        var fullDeferred = $q.defer();

        $q.all(dictPromises)
        .then(function(processedDict) {
            angular.extend(dict, processedDict);

            console.log(dict);

            var stepsChain = $q.when();

            angular.forEach(file.steps, function(step, stepIndex) {

                //run function next in chain (so dict populated)
                stepsChain.then(function() {

                    //new deferred for each step -- or return something in Clotho.run block
                    var stepDeferred = $q.defer();

                    //parse inputs
                    //future - handle promises in inputs
                    //todo - handle Clotho specifics after parse (e.g. enzyme)
                    var inputs = [];
                    angular.forEach(step.input, function(input){

                        //console.log('\n\n\n\ninput:', input);

                        //todo - better handling of empty $parse -- should still show up in computed keys

                        //future - handle 2+ layers deep
                        //for now, go one layer deep
                        if (!angular.isString(input)) {
                            //assume array or object
                            var parsedInput = (angular.isArray(input)) ? [] : {};
                            angular.forEach(input, function(child) {

                                var parsed = $parse(child)(dict) || {};
                                //console.log(child, parsed, dict[child]);

                                //todo - handle object scenario (don't push)
                                //todo - (later) handle defining object keys, e.g. blah[0] or blah.key

                                //for non-dictionary values, or object defined
                                //e.g. blah.value[0]
                                parsed = parsed.value || parsed;

                                parsedInput.push(parsed);
                            });
                            //console.log(parsedInput);
                            inputs.push(parsedInput);
                        } else {
                            parsedInput = $parse(input)(dict) || '';
                            //for non-dictionary values e.g. a boolean or number
                            parsedInput = parsedInput.value || parsedInput;
                            //console.log(parsedInput);
                            inputs.push(parsedInput);
                        }
                    });

                    //console.log('running ' + step.reaction + ' with inputs:', inputs);

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
                    //console.log('result of '+step.reaction+' (#'+stepIndex+') with inputs:', inputs, result);
                    dict[step.output] = {key : step.output, computed : true, stepNum: stepIndex, value : result};
                    stepDeferred.resolve();

                    return stepDeferred.promise;
                });
            });

            return stepsChain;
        })
        .then(function() {

            //console.log('final result: ', dict.final.value, dict);
            fullDeferred.resolve((!!returnDictionary) ? dict : dict.final.value);
        });

        return fullDeferred.promise;
    };

    return {
        process : process
    }


}]);

//todo - returns array, need to maintain object references
Application.Dna.filter('orderByProp', function() {
    return function(obj, prop) {

        //todo - necessary??
        var array = [];
        _.each(obj, function (value) {
            array.push(value);
        });

        return _.sortBy(array, function (obj) {
            return obj.prop
        });
    }
});

Application.Dna.filter('objAddKeys', function() {
    return function(obj) {
        return _.each(obj, function (value, key) {
            _.extend(value, {key : key})
        });
    }
});

Application.Dna.filter('stepCheck', function() {
    return function(input, index) {
        return _.pick(input, function(value, key) {
            return value.computed == false || value.stepNum < index;
        });
    }
});

Application.Dna.directive('constructionDictionaryView', [function() {
    return {
        restrict : 'EA',
        require : 'ngModel',
        scope : {
            dict : '=ngModel',
            editable : '=constructionEditable',
            onChange : '&?constructionOnchange'
        },
        template : '<div style="overflow: scroll">' +
            '<pre class="pre-scrollable" ng-bind="dict | json"></pre>' +
            '<table class="table table-bordered table-condensed constructionDictionary">' +
            '<thead><tr><th>Key</th> <th>Value</th></tr></thead>' +
            '<tbody>' +
                '<tr ng-repeat="item in dict | orderByProp:\'computed\'" ng-class="{\'computed\' : item.computed}">' +
                '<td> <code contenteditable="{{ !!editable }}" ng-model="item.key"></code> </td>' +
                //'<td>{{item}}</td>' +
                '<td contenteditable="{{ !!editable && !item.computed && !item.retrieved}}" ng-model="item.value" style="overflow: hidden; text-overflow: ellipsis;"></td>' +
                '</tr>' +
                '<tr ng-if="!!editable"><td colspan="2"><button class="btn" ng-click="addTerm()">Add new Key</button></td></tr>' +
            '</tbody>' +
            '</table>' +
            '</div>',
        link : function(scope, element, attrs, ngModel) {

            scope.addTerm = function() {
                //custom fields will be added in on process
                //todo -- add ability to clotho parse
                angular.extend(scope.dict, {"" : ""});
            };

            //deep watch
            scope.$watch('dict', function(newval, oldval) {
                console.log(newval, oldval);
                if (!newval || !oldval) return;

                console.log('change to dictionary, running callback');

                angular.isDefined(scope.dict) && scope.onChange({newDict : scope.dict});
            }, true)
        }
    }
}]);

Application.Dna.directive('constructionField', ['$compile', '$filter', function($compile, $filter) {
    return {
        restrict : "EA",
        require: "ngModel",
        scope: {
            fieldName : '@constructionField',
            model : '=ngModel',
            dict : '=?constructionDictionary',
            options : '=constructionOptions',
            required : '@constructionRequired',
            type : '@constructionType'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {

                    var template = angular.element('<div class="control-group" ng-class="{\'error\' : hasInvalidItem}">' +
                        '<label class="control-label">{{ fieldName }}</label>' +
                        '<div class="controls">' +
                        '<constructionFieldInput></constructionFieldInput>' +
                        '</div>' +
                        '</div>');

                    template.find('constructionFieldInput').html(selectTemplate(scope.type));

                    function selectTemplate(type) {
                        switch (angular.lowercase(type)) {
                            case 'value' : {
                                return '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-change="checkInput()" ng-disabled="!options.editable" ng-required="required"/>';
                            }
                            case 'boolean' : {
                                return '<button btn-checkbox class="btn" ng-class="{\'btn-success\' : !!model}"  ng-model="model" ng-disabled="!options.editable"><i ng-class="{\'icon-white icon-ok-circle\' : !!model, \'icon-ban-circle\' : !model}"></i></button>';
                            }
                            case 'array' : {
                                if (angular.isString(scope.model))
                                    scope.model = [scope.model];
                                return '<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-list ng-change="checkInput()" ng-disabled="!options.editable" ng-required="required"/>';
                            }
                            default : {
                                return '<select class="span12" placeholder="{{fieldName}}" ng-model="model" ng-options="key as key for (key,val) in dictFiltered" ng-disabled="!options.editable"  ng-required="required"></select>'
                            }
                        }
                    }

                    element.html($compile(template)(scope));

                },
                post: function postLink(scope, element, attrs, ngModel) {


                    //console.log(scope.options, scope.dict);
                    if (!!scope.dict) {
                        console.log(scope.options.stepIndex, scope.model);
                        _.each(scope.dict, function(val, key) {
                            console.log(key, val.stepNum);
                        });

                        scope.dictFiltered = $filter('stepCheck')(scope.dict, scope.options.stepIndex);
                        console.log('now filtered to:', _.keys(scope.dictFiltered));
                    }


                    scope.required = (scope.required == 'false');

                    scope.hasInvalidItem = false;
                    scope.checkInput = function() {

                        if (scope.type == 'array') {
                            var checkArrayRound = false;
                            angular.forEach(scope.model, function(item) {
                                if (!scope.dictFiltered[item]) {
                                    checkArrayRound = true;
                                }
                            });
                            scope.hasInvalidItem = !!(checkArrayRound);
                        } else {
                            scope.hasInvalidItem = !!scope.dictFiltered[scope.model];
                        }
                    };
                }
            }
        }
    }
}]);

Application.Dna.directive('constructionStep', ['$parse', '$compile', '$http', '$templateCache', '$filter', '$timeout', function($parse, $compile, $http, $templateCache, $filter, $timeout) {
    return {
        restrict : "EA",
        require: "ngModel",
        scope: {
            step : '=ngModel',
            dict : '=constructionDictionary',
            index : '=constructionIndex',
            editable : '@constructionEditable',
            fields : '=constructionStep'
        },
        compile: function compile(tElement, tAttrs, transclude) {
            return {
                pre: function preLink(scope, element, attrs, ngModel) {

                    scope.compileStep = function compileStep() {

                        scope.options = {
                            "editable" : scope.editable,
                            "stepIndex" : scope.index
                        };

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

                            return '<div class="span4" construction-field="'+name+'" construction-type="'+type+'" construction-options="options" construction-dictionary="dict" ng-model="'+model+'" ></div>';
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
                            $http.get('/dna/construction/steps/' + scope.step.reaction + '.html', {cache : $templateCache})
                                .then(function(data) {
                                    //works - great
                                    //console.log('got template for ' + scope.step.reaction + '\n', data.data)
                                    //console.log(data.data);
                                    //console.log(scope);
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
                    };

                    scope.compileStep();

                },
                post: function postLink(scope, element, attrs, ngModel) {

                    //todo -- filter dictionary for step -- apparently need to timeout, may be problematic if run on server and async
                    //scope.dict = $filter('stepCheck')(scope.dict, scope.options.stepIndex);
                    //todo --- watch on index and dictionary in step, pass down trimmed to fields.... OR don't trim dictionary, error class when not in trimmed dictionary OR just check stepNum in ng-class (but set up watch on stepnum, since updates later apparently i.e. after process
                    //alternatively, filter to disable select options
                        // http://jsfiddle.net/alalonde/dZDLg/9/


                    scope.showReaction = true;
                    scope.toggleReaction = function() {
                        scope.showReaction = !scope.showReaction;
                    };

                    scope.$watch('index', function(newind, oldind) {
                        console.log('compiling step ' + newind + ' from ' + oldind);
                        scope.compileStep();
                    })
                }
            }
        }
    }
}]);