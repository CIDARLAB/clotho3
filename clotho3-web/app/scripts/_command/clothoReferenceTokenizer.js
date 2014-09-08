angular.module('clotho.tokenizer')
.directive('clothoReferenceTokenizer',  function ($parse, Clotho, ClientAPI, Debug, clothoTokenCollectionFactory, ClothoSchemas) {

    var Debugger = new Debug('clothoReferenceTokenizer', '#ee9955');

    return {
      restrict: 'E',
      templateUrl: 'views/_command/clothoReferenceTokenizer.html',
      scope: {
        startingTags : '='
      },
      controller: function ($scope, $element, $attrs) {
        this.addToken = function (item) {
          Debugger.log('adding token', item);
          $scope.tokenCollection.addToken(item);
        };

        this.removeToken = function (index, model) {
          Debugger.log('removing token', index);
          $scope.tokenCollection.removeToken(index);
        };

        this.handleBackout = function () {
          var activeTokenIndex = $scope.tokenCollection.whichActive();
          if (activeTokenIndex >= 0) {

          } else {
            $scope.tokenCollection.setLastActive()
          }
        };

        this.toggleTokenActive = function (index, event) {
          event.preventDefault();
          $scope.tokenCollection.toggleActive(index);
        };

        this.handleSelect = function (item, query) {
          this.addToken(item);
        };

        this.disambiguate = function () {
          console.log($scope.tokenCollection.tokens);

          //todo - construct run if function, or show view if not

          //get first token
          var isFunction = false;

          ClientAPI.say({
            from: 'client',
            channel: 'run',
            class: 'info',
            text: 'disambiguating...'
          });

          if (isFunction) {
            //todo
            Clotho.run().then(function (response) {
              ClientAPI.say({
                from: 'server',
                channel: 'submit',
                class: 'success',
                text: response
              });
              $scope.query = '';
            });
          } else {

          }
        };

        this.focusInput = function () {
          $element.find('[clotho-reference-autocomplete]').focus();
        };
      },
      controllerAs: 'tokenCtrl',
      link: function (scope, element, attrs) {

        /* Tokens */

        scope.tokenCollection = new clothoTokenCollectionFactory(scope.startingTags);

        scope.$watchCollection('tokenCollection.tokens', function () {
          // todo - check type, should
          // (1) filter next token if function OR
          // (2) prevent more tokens if sharable

          //todo - allow primitives if schema not specified (or wants string)

          scope.currentFilter = function (item) {
            return ClothoSchemas.isFunction(item);
          };
        });
      }
    };
});

