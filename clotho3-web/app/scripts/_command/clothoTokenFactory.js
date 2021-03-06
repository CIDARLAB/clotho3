angular.module('clotho.tokenizer')
/**
 * @name clothoTokenFactory
 *
 * @description
 * ClothoTokens are essentially wrappers for clotho sharables, or strings.
 * They expect the fields minimally of name, id, schema to be a sharable, or just a string for other keywords
 */
  .factory('clothoTokenFactory', function ($q, Clotho) {

    //pass UUID to make object, or just pass value as string
    function ClothoToken(sharable) {
      var self = this;
      self.model = sharable;
      self.activeState = false;

      if (this.isSharable()) {
        self.fullSharablePromise = Clotho.get(self.model.id, {mute: true}).then(function (data) {
          self.fullSharable = data;
          return data;
        });
      } else {
        self.fullSharablePromise = $q.when(null);
      }
    }

    ClothoToken.prototype.readable = function () {
      return this.model.name || this.model;
    };

    ClothoToken.prototype.id = function () {
      return this.model.id || null;
    };

    ClothoToken.prototype.isAmbiguous = function () {
      return angular.isArray(this.model);
    };

    ClothoToken.prototype.isSharable = function () {
      return !this.isAmbiguous() && angular.isDefined(this.model.id);
    };

    //getter setter
    ClothoToken.prototype.active = function (state) {
      if (angular.isDefined(state)) {
        this.activeState = !!state;
      } else {
        return this.activeState === true;
      }
    };

    return ClothoToken;
  });
