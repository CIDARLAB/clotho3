'use strict';

/**
 * @name unique
 * @example {{items | unique:attribute | json}}
 * Filters out all duplicate items from an array by checking the specified key
 * @param [key] {string} the name of the attribute of each object to compare for uniqueness
 if the key is empty, the entire object will be compared
 if the key === false then no filtering will be performed
 * @return {array}
 */
Application.Interface.filter('unique', ['$parse', function ($parse) {

    return function (items, filterOn) {

        if (filterOn === false) {
            return items;
        }

        if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
            var hashCheck = {}, newItems = [],
                get = angular.isString(filterOn) ? $parse(filterOn) : function (item) { return item; };

            var extractValueToCompare = function (item) {
                return angular.isObject(item) ? get(item) : item;
            };

            angular.forEach(items, function (item) {
                var valueToCheck, isDuplicate = false;

                for (var i = 0; i < newItems.length; i++) {
                    if (angular.equals(extractValueToCompare(newItems[i]), extractValueToCompare(item))) {
                        isDuplicate = true;
                        break;
                    }
                }
                if (!isDuplicate) {
                    newItems.push(item);
                }

            });
            items = newItems;
        }
        return items;
    };
}]);

/**
 * @name highlight
 * @example

 <pre class="prettyprint">
 &lt;label&gt;&lt;input type=&quot;checkbox&quot; ng-model=&quot;caseSensitive&quot;&gt; Case Sensitive?&lt;/label&gt;
 &lt;input placeholder=&quot;Enter some text to highlight&quot; value=&quot;you&quot; ng-model=&quot;highlightText&quot;&gt;
 &lt;p ng-bind-html-unsafe=&quot;&#x27;Hello there, how are you today? I\'m fine thank you.&#x27; | highlight:highlightText:caseSensitive&quot;&gt;&lt;/p&gt;

 &lt;style&gt;
 .ui-match { background: yellow; }
 &lt;/style&gt;
 </pre>

 * @param text {string} haystack to search through
 * @param search {string} needle to search for
 * @param [caseSensitive] {boolean} optional boolean to use case-sensitive searching
 */
Application.Interface.filter('highlight', function () {
    return function (text, search, caseSensitive) {
        if (search || angular.isNumber(search)) {
            text = text.toString();
            search = search.toString();
            if (caseSensitive) {
                return text.split(search).join('<span class="ui-match">' + search + '</span>');
            } else {
                return text.replace(new RegExp(search, 'gi'), '<span class="ui-match">$&</span>');
            }
        } else {
            return text;
        }
    };
});


/**
 * A replacement utility for internationalization very similar to sprintf.
 *
 * @param replace {mixed} The tokens to replace depends on type
 *  string: all instances of $0 will be replaced
 *  array: each instance of $0, $1, $2 etc. will be placed with each array item in corresponding order
 *  object: all attributes will be iterated through, with :key being replaced with its corresponding value
 * @return string
 *
 * @example: 'Hello :name, how are you :day'.format({ name:'John', day:'Today' })
 * @example: 'Records $0 to $1 out of $2 total'.format(['10', '20', '3000'])
 * @example: '$0 agrees to all mentions $0 makes in the event that $0 hits a tree while $0 is driving drunk'.format('Bob')
 */
Application.Interface.filter('format', function(){
    return function(value, replace) {
        var target = value;
        if (angular.isString(target) && replace !== undefined) {
            if (!angular.isArray(replace) && !angular.isObject(replace)) {
                replace = [replace];
            }
            if (angular.isArray(replace)) {
                var rlen = replace.length;
                var rfx = function (str, i) {
                    i = parseInt(i, 10);
                    return (i>=0 && i<rlen) ? replace[i] : str;
                };
                target = target.replace(/\$([0-9]+)/g, rfx);
            }
            else {
                angular.forEach(replace, function(value, key){
                    target = target.split(':'+key).join(value);
                });
            }
        }
        return target;
    };
});
