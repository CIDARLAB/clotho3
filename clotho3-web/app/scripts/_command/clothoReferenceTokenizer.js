angular.module('clotho.tokenizer')
.directive('clothoReferenceTokenizer',  function ($parse, Clotho, ClientAPI, Debug, clothoTokenCollectionFactory, ClothoSchemas, $timeout) {

    var Debugger = new Debug('clothoReferenceTokenizer', '#ee9955');

    return {
      restrict: 'E',
      templateUrl: 'views/_command/clothoReferenceTokenizer.html',
      scope: {
        startingTags : '='
      },
      controller: function ($scope, $element, $attrs) {

        function resetQuery () {
          $scope.query = '';
        }

        function resetTokens() {
          $scope.tokenCollection.removeAll();
        }

        this.addToken = function (item, isPrimitive) {
          Debugger.log('adding token', item);
          $scope.tokenCollection.addToken(item);
        };

        this.removeToken = function (index, model) {
          Debugger.log('removing token', index);
          $scope.tokenCollection.removeToken(index);
        };

        this.handleBackout = function (evt) {
          if ($scope.tokenCollection.isLastActive()) {
            $scope.tokenCollection.clearLast();
          } else {
            //hack - timeout to allow tokenCollection event handlers
            $timeout(function() {
              $scope.tokenCollection.setLastActive();
            });
          }
        };

        this.onTokenActive = function (index, event, token) {

        };

        this.handleSelect = function (item, query) {
          //select will not get called unless there is a dropdown, check enter callback
          this.addToken(item);
        };

        this.handleEnter = function (event, query) {
          event.preventDefault();
          if ($scope.allowPrimitive === true && query.length) {
            this.addToken(query);
          } else {
            this.disambiguate(event)
          }
        };

        this.disambiguate = function (event) {

          var tokens = $scope.tokenCollection.tokens;
          if (tokens.length) {
            var tokenNames = [],
              firstToken = tokens[0];
            angular.forEach(tokens, function (token) {
              tokenNames.push(token.readable());
            });

            ClientAPI.say({
              from: 'client',
              channel: 'run',
              class: 'info',
              tokens: angular.copy(tokens), //prevent reference updates
              text: 'disambiguating: ' + tokenNames.join(' ')
            });

            //parsing... run function, show View, edit sharable / schema

            if (ClothoSchemas.isFunction(firstToken.model)) {

              var args = angular.map(tokens.slice(1), function (token) {
                return token.isSharable() ? token.model.id : token.model;
              });

              Clotho.run(firstToken.model.id, args).then(function (response) {
                ClientAPI.say({
                  from: 'server',
                  channel: 'submit',
                  class: 'success',
                  text: response
                });
              });
            }
            //check if view
            else if (ClothoSchemas.isView(tokens[0])) {
              console.log('need to handle view from command bar');
            }
            //not a function
            else {
              Clotho.edit(tokens[0].id);
            }

            resetQuery();
            resetTokens();
          }
          //no tokens
          else {
            //do nothing
          }
        };

        this.focusInput = function () {
          $element.find('[clotho-reference-autocomplete]').focus();
        };
      },
      controllerAs: 'tokenCtrl',
      link: function (scope, element, attrs) {

        function clearFilter () {
          scope.currentFilter = null;
          scope.currentPlaceholder = null;
        }

        function setFilter (filter) {
          scope.currentFilter = filter;
        }

        function setPlaceholder (placeholder) {
          scope.currentPlaceholder = placeholder;
        }

        function allowPrimitive (val) {
          scope.allowPrimitive = angular.isDefined(val) ? !!val : true;
        }

        function disallowPrimitive () {
          scope.allowPrimitive = false;
        }

        function setHideFilter () {
          setFilter(function () {
            return false;
          });
        }

        function blockInput () {
          setHideFilter();
          scope.blockInput = true;
        }

        function unblockInput () {
          scope.blockInput = false;
          clearFilter();
        }

        /* Tokens */

        scope.tokenCollection = new clothoTokenCollectionFactory(scope.startingTags);

        scope.$watchCollection('tokenCollection.tokens', function (newval, oldval) {

          //if we have tokens
          if (newval && newval.length > 0) {

            var firstToken = newval[0], //ClothoToken class
              numTokens = newval.length;

            //if first token is a function, inspect arguments
            //set filter for current token based on function args + numTokens
            if (ClothoSchemas.isFunction(firstToken.model)) {

              firstToken.fullSharablePromise.then(function (fullFunction) {

                //check that args are defined
                if (fullFunction.args) {

                  //if we have more arguments left
                  if (fullFunction.args.length > (numTokens - 1)) {
                    var curArg = fullFunction.args[numTokens - 1];

                    //todo - handle arg.type undefined

                    var isPrimitive = ClothoSchemas.isPrimitiveField(curArg.type);

                    if (isPrimitive) {
                      setHideFilter();
                      setPlaceholder('Enter ' + curArg.type);
                      allowPrimitive();
                    } else {
                      setFilter(function (item) {
                        return ClothoSchemas.isInstanceOfSchema(item, curArg.type);
                      });
                      setPlaceholder(curArg.type);
                      disallowPrimitive();
                    }
                  }
                  //arguments are full
                  else {
                    blockInput();
                    setPlaceholder('All arguments defined');
                    disallowPrimitive();
                  }
                }
                //no function args defined... so allow any input
                else {
                  unblockInput();
                  allowPrimitive();
                }
              });
            }
            //first token not a function
            else {
              blockInput();
              setPlaceholder('Does not accept arguments');
            }
          }
          //no tokens are defined... no filter
          else {
            unblockInput();
            disallowPrimitive();
          }
        });
      }
    };
});

