angular.module('clotho.tokenizer')
/**
 * @name clothoTokenCollectionFactory
 *
 * @description
 * Object to handle a collection of ClothoTokens
 *
 * todo - active state storage in token, not collection
 */
	.factory('clothoTokenCollectionFactory', function (clothoTokenFactory) {

    //want this easily accessible globally
    //get token at index, otherwise last one
		var getToken = angular.noop;

    function ClothoTokenCollection (startingTokens) {

      var self = this;
			this.tokens = [];

			if (angular.isArray(startingTokens) || angular.isObject(startingTokens)) {
				angular.forEach(startingTokens, function (token) {
					self.addToken(token);
				});
			} else {
        //wrong format
      }

      getToken = function getToken (index) {
        if (angular.isDefined(index)) {
          return self.tokens[index];
        } else {
          return self.tokens[self.tokens.length - 1];
        }
      }
		}

		//add a token, pass arguments through to clothoTokenFactory
		ClothoTokenCollection.prototype.addToken = function (sharable) {
			this.tokens.push(new clothoTokenFactory(sharable));
		};

    ClothoTokenCollection.prototype.hasLength = function () {
      return this.tokens.length > 0;
    }

		ClothoTokenCollection.prototype.inRange = function (index) {
			return index > -1 && index < this.tokens.length;
		};

		//get token at given index
		ClothoTokenCollection.prototype.getToken = function (index) {
			return getToken(index);
		};

		//return index of token
		ClothoTokenCollection.prototype.indexOf = function (token) {
			return this.tokens.indexOf(token);
		};

		//remove token at given index, return it if removed, otherwise false
		ClothoTokenCollection.prototype.removeToken = function (index) {
			if (this.inRange(index)) {
        //todo - check memory leak (shouldn't have reference but...)
				return this.tokens.splice(index, 1);
			} else {
				return false;
			}
		};

		//remove all tokens
		ClothoTokenCollection.prototype.removeAll = function () {
			this.tokens.length = 0;
		};

		//remove all active tokens. return array of removed tokens
		ClothoTokenCollection.prototype.removeActiveTokens = function () {
      var removed = [];
      for (var i = 0; i < this.tokens.length; i++) {
        if (getToken(i).active()) {
          var toReturn = this.removeToken(i);
          removed.push(toReturn);
        }
      }
      return removed;
		};

		//set token at given index to be active
		ClothoTokenCollection.prototype.setActive = function (index) {
			if (this.inRange(index)) {
				getToken(index).active(true);
        return true;
			} else {
				return false;
			}
		};

		//unset active tokens
		ClothoTokenCollection.prototype.unsetActive = function () {
			angular.forEach(this.tokens, function (token) {
        token.active(false);
      });
		};

		//check if token is active, at index when passed, otherwise if any is active
		ClothoTokenCollection.prototype.isActive = function (index) {
			if (angular.isDefined(index)) {
				return getToken(index).active();
			} else {
				for (var i = 0; i < this.tokens.length; i++) {
          if (getToken(i).active()) {
            return true;
          }
        }
        return false;
			}
		};

    //set token at last position to be active
    ClothoTokenCollection.prototype.setLastActive = function () {
      return this.hasLength() && getToken().active(true);
    };

		ClothoTokenCollection.prototype.isLastActive = function () {
			return this.hasLength() && getToken().active();
		};

    ClothoTokenCollection.prototype.clearLast = function () {
      return this.removeToken(this.tokens.length - 1);
    };

		return ClothoTokenCollection;
	});
