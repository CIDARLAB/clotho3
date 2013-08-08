/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import lombok.Getter;

/**
 *
 * @author spaige
 */
    public final class Argument{
        @Getter
        private final String name;
        @Getter
        private final Class type;
        public Argument(String name, Class type){
            this.name = name;
            this.type = type;
        }
    }